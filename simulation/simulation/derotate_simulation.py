import os
import simulation
import numpy as np
import pandas as pd
from simulation.parsers.parse_trace import parse_trace
from simulation.units import gadget_time_units
from simulation.util import get_sim_name, get_pivot, get_quat, get_omega_mb
import pynbody
from astropy.table import Table
import matplotlib.pyplot as plt
from pyquaternion import Quaternion
import tqdm
import quaternion
import argparse
import gc


def rotate_vec(vec, quat):
    """Rotate a numpy array of 3-vectors `vec` given a quaternion `quat`"""
    v = np.hstack([np.zeros((len(vec), 1)), vec])
    vq = quaternion.as_quat_array(v)
    new_vec = quat * vq * quat.conj()

    # remove w component
    return quaternion.as_float_array(new_vec)[:, 1:]


def rotate_simarray(vec, quat):
    """Rotate a SimArray array of 3-vectors `vec` given a quaternion `quat`"""
    new_vec = rotate_vec(vec, quat)
    return pynbody.array.SimArray(new_vec, vec.units)


def rotate_on_orbit_plane(pos, vel):
    # Rotation of -90deg around x axis
    orbit_plane_quat = 1/np.sqrt(2) * np.quaternion(1, -1, 0, 0)
    new_pos = rotate_simarray(pos, orbit_plane_quat)
    new_vel = rotate_simarray(vel, orbit_plane_quat)
    return new_pos, new_vel


def get_all_keys(snap):
    """return all the (non derived) keys for all the families"""
    ak = set()
    for fam in snap.families():
        ak = ak.union(snap[fam].loadable_keys()).union(snap[fam].keys()).union(snap[fam].family_keys())
    ak = [k for k in ak if not k in ['pos', 'vel', 'acce', "x", "y", "z", "vx", "vy", "vz", 'acce_x', 'acce_y', 'acce_z']]
    ak.sort()
    return ak

def rotate_snap(input_snap, quat, omega_mb, pivot, offset=None, on_orbit_plane=False):
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

    new_pos, new_vel = derotate_pos_and_vel(f['pos'], f['vel'], quat, omega_mb, pivot)

    if on_orbit_plane:
        print("Rotating on the plane of the orbit...")
        new_pos, new_vel = rotate_on_orbit_plane(new_pos, new_vel)

    del s['pos']
    s['pos'] = new_pos.astype(f['pos'].dtype)

    del s['vel']
    s['vel'] = new_vel.astype(f['vel'].dtype)

    # print(f['pos'])
    # print(s['pos'])
    # print(f.properties)

    s.properties = f.properties.copy()
    s.properties['z'] = f.properties['z']
    s.properties['boxsize'] = 60 * pynbody.units.kpc  ## FIXME why not copying boxsize of f
    # print(s.properties)
    return s


def derotate_pos_and_vel(pos, vel, quat, omega_mb, pivot, offset=None):

    # pos
    new_pos = rotate_vec(pos - pivot, quat)
    if offset is not None:
        new_pos += offset
    new_pos = pynbody.array.SimArray(new_pos, pos.units)

    # vel
    new_vel = rotate_vec(vel + np.cross(omega_mb, pos - pivot), quat)
    new_vel = pynbody.array.SimArray(new_vel, vel.units)
    return new_pos, new_vel


def write_rotated_snap(s, filename):
    s.write(fmt=pynbody.snapshot.gadget.GadgetSnap, filename=filename, use_time=True)


def get_quaternions(trace):
    quat = np.vstack([trace.quat_w, trace.quat_x, trace.quat_y, trace.quat_z]).T
    return quaternion.as_quat_array(quat)


def derotate_simulation(sim_path, new_path, snap_indexes=slice(None, None, None), on_orbit_plane=False):
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

    for i, (snap, q, om_mb) in enumerate(tqdm.tqdm(zip(sim, loc_quat, loc_omega_mb), total=len(sim))):
        # print(snap)
        # print(quat)
        s_rot = rotate_snap(snap, q, om_mb, pivot, offset[i], on_orbit_plane)
        # This requires a change in source code of pynbody since it seems not possible
        # to set the dtype of the initial arrays as created by new()
        assert s_rot['mass'].dtype == np.float32
        assert s_rot['pos'].dtype == np.float32
        assert s_rot['vel'].dtype == np.float32

        write_rotated_snap(s_rot, os.path.join(new_path, os.path.basename(snap._filename)))
        del snap
        sim.snap_list[i] = None  # destroying references to the snap and the list
        if i % 10 == 0:
            gc.collect()

def parse_args(cli=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(dest='sim_path', help="Path to the simulation snapshots")
    parser.add_argument("-o", '--outpath', help="Folder of the derotated snapshots", default=None)
    parser.add_argument('--on-orbit-plane', help='Rotate snaps to be viewed from the orbit plane', action='store_true')
    parser.add_argument('--start', help="First snap index", default=None, type=int)
    parser.add_argument('--stop', help="Last snap index", default=None, type=int)
    parser.add_argument('-n', help="Take one every n snapshots", default=None, type=int)
    args = parser.parse_args(cli)
    return args


def main(cli=None):
    args = parse_args(cli)
    sim_name = get_sim_name(args.sim_path)
    if args.outpath is None:
        new_path = sim_name + "_derot"
    else:
        new_path = args.outpath
    derotate_simulation(args.sim_path, new_path, slice(args.start, args.stop, args.n), on_orbit_plane=args.on_orbit_plane)


if __name__ == '__main__':
    main()