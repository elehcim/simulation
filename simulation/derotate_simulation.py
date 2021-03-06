import os
import simulation
import numpy as np
import pandas as pd
from simulation.parsers.parse_trace import parse_trace
from simulation.units import gadget_time_units
from simulation.util import get_sim_name, get_pivot, get_quat, get_omega_mb, get_quat_omega_pivot, get_initial_rotation
import pynbody
from astropy.table import Table
import matplotlib.pyplot as plt
import tqdm
import quaternion
import argparse
import gc

ORBIT_PLANE_QUAT = 1/np.sqrt(2) * np.quaternion(1, -1, 0, 0)


class Derotator:
    def __init__(self, sim, sim_name=None):
        """Create sliced tables of quaternions, omega. Get pivot too
        They are sliced so that they can be indexed with an enumerate list on `sim`"""
        my_sim_name = sim_name or get_sim_name(sim.sim_id)
        # logger.info(f"Initialize derotator from {my_sim_name}...")

        slicer = sim._snap_indexes
        quat_arr, omega_mb_arr, pivot = get_quat_omega_pivot(my_sim_name)

        self.my_sim_name = my_sim_name
        self.quat_arr = quat_arr[slicer]
        self.omega_mb_arr = omega_mb_arr[slicer]
        self.omega_mb_0 = self.omega_mb_arr[0]
        self.pivot = pivot
        # print(self.quat_arr, self.omega_mb_arr, self.pivot)

        # try:
        #     quat_vp0 = get_initial_rotation(initial_rotation_simname)
        #     logger.info('Apply initial rotation...')

        #     print(quat_vp0)
        # except Exception:
        #     logger.warning('Cannot get initial rotation.')

    def derotate_snap(self, snap, idx_snap, _initial_rotation=False):
        """I dont like it, maybe a dict like structure would be better, but for now I leave it like that.
        The assumpion here is to have the index synchronized with the slicer (snap_indexes) of the sim
        """
        quat = np.quaternion(*self.quat_arr[idx_snap])
        omega_mb = self.omega_mb_arr[idx_snap]# - self.omega_mb_0

        return derotate_pos_and_vel(snap['pos'], snap['vel'], quat, omega_mb, self.pivot)

    def derotate_sim(self, sim, on_orbit_plane, _initial_rotation=False):
        assert len(self.quat_arr) == len(sim)
        assert len(self.omega_mb_arr) == len(sim)
        for i, snap in enumerate(tqdm.tqdm(sim)):
            quat = np.quaternion(*self.quat_arr[i, :])
            omega_mb = self.omega_mb_arr[i, :]
            # logger.info("Derotating...")
            # logger.debug("quat:     {}".format(quat))
            # logger.debug("omega_mb: {}".format(omega_mb))
            # logger.debug("pivot:    {}".format(self.pivot))
            snap['pos'], snap['vel'] = derotate_pos_and_vel(snap['pos'], snap['vel'], quat, omega_mb, self.pivot)


        # This is important in order to compare angular momentum from Moria to MovingBox
        # It is important to do this before the on_orbit_plane rotation
        if _initial_rotation:
            # logger.info('Apply initial rotation...')
            quat_vp0 = get_initial_rotation(self.my_sim_name)
            print(quat_vp0)
            snap['pos'] = rotate_vec(snap['pos'], quat_vp0)
            snap['vel'] = rotate_vec(snap['vel'], quat_vp0)


        if on_orbit_plane:
            # logger.debug("Rotating on the plane of the orbit...")
            snap['pos'], snap['vel'] = rotate_on_orbit_plane(snap['pos'], snap['vel'])



def rotate_vec(vec, quat):
    """Rotate a numpy array of 3-vectors `vec` given a quaternion `quat`
    Since I need to get the position (X*) in the absolute reference frame from the position in the box (x*),
    the formula is X* = q  x*  q^-1 (the inverse of eq. (1) in Nichols 2015).
    This is because their quaternion is to go from the fixed (a.k.a. "inertial", "host") galaxy frame to the periodic box.
    My quaternion is to rotate a vector in the box to the one in the fixed frame.

    Note that this returns a vector with shape vec.shape + quat.shape, i.e. all the vec rows are rotated by all the quat rows.
    If you want to rotate each vec row for each corresponding quat row (needless to say that len(vec) == len(quat)), you want
    to use the rotate_vec_by_quat_array function. FIXME, there should be a way to do it in one shot.

    """

    assert quat.dtype == np.quaternion, "Please provide a np.quaternion as input to rotate_vec"
    return quaternion.rotate_vectors(quat, vec)

def rotate_vec_by_quat_array(vec, quat):
    """Rotate a numpy array of 3-vectors `vec` given a quaternion `quat`
    Since I need to get the position (X*) in the absolute reference frame from the position in the box (x*),
    the formula is X* = q  x*  q^-1 (the inverse of eq. (1) in Nichols 2015).
    This is because their quaternion is to go from the fixed (a.k.a. "inertial", "host") galaxy frame to the periodic box.
    My quaternion is to rotate a vector in the box to the one in the fixed frame.
    """

    assert quat.dtype == np.quaternion, "Please provide a np.quaternion as input to rotate_vec"
    v = np.hstack([np.zeros((len(vec), 1)), vec])
    vq = quaternion.as_quat_array(v)
    new_vec = quat * vq * quat.conj()

    # In test_quaternion.py I verified that the above rotation quat * vq * quat.conj() is equivalent to quaternion.rotate_vectors(quat, vec)

    # remove w component
    return quaternion.as_float_array(new_vec)[:, 1:]


def rotate_simarray(vec, quat):
    """Rotate a SimArray array of 3-vectors `vec` given a quaternion `quat`"""
    new_vec = rotate_vec(vec, quat)
    return pynbody.array.SimArray(new_vec, vec.units)

def rotate_simarray_on_orbit_plane(vec, orbit_plane_quat=ORBIT_PLANE_QUAT):
    # Rotation of -90deg around x axis
    new_vec = rotate_simarray(vec, orbit_plane_quat)
    return new_vec

def rotate_vec_on_orbit_plane(vec, orbit_plane_quat=ORBIT_PLANE_QUAT):
    # Rotation of -90deg around x axis
    new_vec = rotate_vec(vec, orbit_plane_quat)
    return new_vec

def rotate_on_orbit_plane(pos, vel, orbit_plane_quat=ORBIT_PLANE_QUAT):
    # Rotation of -90deg around x axis
    new_pos = rotate_simarray(pos, orbit_plane_quat)
    new_vel = rotate_simarray(vel, orbit_plane_quat)
    return new_pos, new_vel


def get_all_keys(snap):
    """return all the (non derived) keys for all the families"""
    ak = set()
    for fam in snap.families():
        ak = ak.union(snap[fam].loadable_keys()).union(snap[fam].keys()).union(snap[fam].family_keys())
    # These are the common keys for all the families.
    ak = [k for k in ak if not k in ['mass', 'pos', 'vel', 'acce', "x", "y", "z", "vx", "vy", "vz", 'acce_x', 'acce_y', 'acce_z']]
    ak.sort()
    return ak

def rotate_snap(input_snap, quat, omega_mb, pivot, on_orbit_plane=False, offset=None):
    f = input_snap
    s = pynbody.snapshot.new(dm=len(f.dm), gas=len(f.gas), star=len(f.star), order='gas,dm,star')
    # print(get_all_keys(f))

    # print(get_all_keys(f.gas))
    for k in get_all_keys(f.gas):
        s.gas[k] = f.gas[k]

    # print(get_all_keys(f.dm))
    for k in get_all_keys(f.dm):
        s.dm[k] = f.dm[k]

    # print(get_all_keys(f.s))
    for k in get_all_keys(f.s):
        s.s[k] = f.s[k]

    new_pos, new_vel = derotate_pos_and_vel(f['pos'], f['vel'], quat, omega_mb, pivot, offset)

    if on_orbit_plane:
        print("Rotating on the plane of the orbit...")
        new_pos, new_vel = rotate_on_orbit_plane(new_pos, new_vel)

    # fix dtypes. The problem seems that in pynbody.new, dtype is assigned to be np.float64 and can't be changed without casting
    del s['pos']
    s['pos'] = new_pos.astype(f['pos'].dtype)

    del s['vel']
    s['vel'] = new_vel.astype(f['vel'].dtype)

    del s['mass']
    s['mass'] = f['mass']

    # print(f['pos'])
    # print(s['pos'])
    # print(f.properties)

    s.properties = f.properties.copy()
    s.properties['z'] = f.properties['z']
    s.properties['boxsize'] = 10000 * pynbody.units.kpc  ## FIXME why not copying boxsize of f
    s.header = f.header
    # print(s.properties)
    return s


def derotate_pos_and_vel(pos, vel, quat, omega_mb, pivot, offset=None):
    """Implementing Eq. 2 in Nichols 2015, the inverse of it actually with the inverse quaternion."""
    # pos
    new_pos = rotate_vec(pos - pivot, quat)
    if offset is not None:
        new_pos += offset
    new_pos = pynbody.array.SimArray(new_pos, pos.units)

    # vel
    new_vel = rotate_vec(vel + np.cross(omega_mb, pos - pivot), quat)
    new_vel = pynbody.array.SimArray(new_vel, vel.units)
    return new_pos, new_vel


def write_rotated_snap(s, filename, out_mat=False):
    if out_mat:
        from simulation.convert_snapshot_to_mat import convert_to_mat_all_info
        convert_to_mat_all_info(s, density_threshold=0.0, outfile_name=filename+".mat")
    else:
        s.write(fmt=pynbody.snapshot.gadget.GadgetSnap, filename=filename, use_time=True)


def get_quaternions(trace):
    quat = np.vstack([trace.quat_w, trace.quat_x, trace.quat_y, trace.quat_z]).T
    return quaternion.as_quat_array(quat)


def derotate_simulation(sim_path, new_path, snap_indexes=slice(None), on_orbit_plane=False, out_mat=False):
    sim = simulation.Simulation(sim_path, snap_indexes=snap_indexes)

    sim_name = get_sim_name(sim_path)
    pivot = get_pivot(sim_name)
    print("Pivot:", pivot)

    locations = np.digitize(sim.times_header, sim.trace.t, right=True)

    quat = get_quat(sim_name)
    omega_mb = get_omega_mb(sim_name)
    if quat is None or omega_mb is None:
        raise RuntimeError('Cannot find quaternion. First save it using save_quat.py')
        # quaternions = get_quaternions(sim.trace)
        # loc_quat = quaternions[locations]
    else:
        loc_quat = quaternion.as_quat_array(quat)[snap_indexes]
        loc_omega_mb = omega_mb[snap_indexes]

    # print(loc_quat)
    print(loc_quat.shape)
    offset = (np.vstack([sim.trace.x, sim.trace.y, sim.trace.z]).T)[locations]

    os.makedirs(new_path, exist_ok=True)
    assert len(loc_quat) == len(sim)
    for i, (snap, q, om_mb) in enumerate(tqdm.tqdm(zip(sim, loc_quat, loc_omega_mb), total=len(sim))):
        # print(snap)
        # print(q)
        # print(pivot)
        s_rot = rotate_snap(snap, q, om_mb, pivot, on_orbit_plane, offset[i])
        # This requires a change in source code of pynbody since it seems not possible
        # to set the dtype of the initial arrays as created by new()
        assert s_rot['mass'].dtype == np.float32
        assert s_rot['pos'].dtype == np.float32
        assert s_rot['vel'].dtype == np.float32

        out_name = os.path.join(new_path, os.path.basename(snap._filename))
        write_rotated_snap(s_rot, out_name, out_mat)
        del snap
        sim.snap_list[i] = None  # destroying references to the snap and the list
        if i % 10 == 0:
            gc.collect()

def parse_args(cli=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(dest='sim_path', help="Path to the simulation snapshots")
    parser.add_argument("-o", '--outpath', help="Folder of the derotated snapshots", default=None)
    parser.add_argument('--on-orbit-plane', help='Rotate snaps to be viewed from the orbit plane', action='store_true')
    parser.add_argument('--out-mat', help='Output matfile', action='store_true')
    parser.add_argument('--start', help="First snap index", default=None, type=int)
    parser.add_argument('--stop', help="Last snap index", default=None, type=int)
    parser.add_argument('-n', help="Take one every n snapshots", default=None, type=int)
    parser.add_argument('--snap_idxs', nargs='+', help="Take these snapshots", default=None, type=int)
    args = parser.parse_args(cli)
    return args


def main(cli=None):
    args = parse_args(cli)
    sim_name = get_sim_name(args.sim_path)
    if args.outpath is None:
        new_path = sim_name + "_derot"
    else:
        new_path = args.outpath
    if args.snap_idxs is not None:
        print("Getting snapshots: ", args.snap_idxs)
        derotate_simulation(args.sim_path, new_path, args.snap_idxs,
                            on_orbit_plane=args.on_orbit_plane, out_mat=args.out_mat)
    else:
        derotate_simulation(args.sim_path, new_path, slice(args.start, args.stop, args.n),
                            on_orbit_plane=args.on_orbit_plane, out_mat=args.out_mat)


if __name__ == '__main__':
    main()