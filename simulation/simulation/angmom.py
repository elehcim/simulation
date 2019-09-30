import pynbody
from pynbody.analysis.angmom import calc_faceon_matrix, calc_sideon_matrix

from astropy import units as u
from simulation.util import get_initial_rotation, get_quat_omega_pivot, get_sim_name, setup_logger
from simulation.derotate_simulation import (rotate_vec, derotate_pos_and_vel, derotate_simulation,
                                           rotate_vec_on_orbit_plane, rotate_on_orbit_plane)
from simulation.units import gadget_angmom_units
from astropy.table import Table
from collections import defaultdict

import gc
import tqdm
import numpy as np
import quaternion

logger = setup_logger('angmom', logger_level='INFO')


def sideon(h, vec_to_xform=calc_sideon_matrix,
             cen=None, vcen=None, move_all=True, **kwargs):
    """
    Reposition and rotate the simulation containing the halo h to see
    h's disk edge on.

    Given a simulation and a subview of that simulation (probably the
    halo of interest), this routine centers the simulation and rotates
    it so that the disk lies in the x-z plane. This gives a side-on
    view for SPH images, for instance.

    """

    from pynbody.analysis.angmom import ang_mom_vec
    from pynbody.analysis import halo
    from pynbody import filt, units, transformation

    if move_all:
        top = h.ancestor
    else:
        top = h

    # Top is the top-level view of the simulation, which will be
    # transformed

    if cen is None:
        # logger.info("Finding halo center...")
        # or h['pos'][h['phi'].argmin()]
        cen = halo.center(h, retcen=True, **kwargs)
        # logger.info("... cen=%s" % cen)

    tx = transformation.inverse_translate(top, cen)

    if vcen is None:
        vcen = halo.vel_center(h, retcen=True)

    tx = transformation.inverse_v_translate(tx, vcen)

    # logger.info("Calculating angular momentum vector...")
    am = ang_mom_vec(h)
    # print(" before rotation L:", am)
    trans = vec_to_xform(am)


    # logger.info("Transforming simulation...")

    tx = transformation.transform(tx, trans)
    # print(" after rotation  L:", ang_mom_vec(h))

    # logger.info("...done!")

    return tx


def faceon(h, **kwargs):
    return sideon(h, vec_to_xform=calc_faceon_matrix, **kwargs)

class Derotator:
    def __init__(self, sim, sim_name=None):
        """Create sliced tables of quaternions, omega. Get pivot too
        They are sliced so that they can be indexed with an enumerate list on `sim`"""
        my_sim_name = sim_name or get_sim_name(sim.sim_id)
        logger.info(f"Initialize derotator from {my_sim_name}...")

        slicer = sim._snap_indexes
        quat_arr, omega_mb_arr, pivot = get_quat_omega_pivot(my_sim_name)

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

    # @static
    # def apply_initial_rotation(cls, initial_rotation_simname):
    #     L = rotate_vec(L, quat_vp0)
    #     Lg = rotate_vec(Lg, quat_vp0)


    def derotate_snap(self, snap, idx_snap, _initial_rotation=False):
        """I dont like it, maybe a dict like structure would be better, but for now I leave it like that.
        The assumpion here is to have the index synchronized with the slicer (snap_indexes) of the sim
        """
        quat = np.quaternion(*self.quat_arr[idx_snap])
        omega_mb = self.omega_mb_arr[idx_snap]# - self.omega_mb_0

        return derotate_pos_and_vel(snap['pos'], snap['vel'], quat, omega_mb, self.pivot)

    def derotate_sim(self, sim, on_orbit_plane):
        assert len(self.quat_arr) == len(sim)
        assert len(self.omega_mb_arr) == len(sim)
        for i, snap in enumerate(tqdm.tqdm(sim)):
            quat = np.quaternion(*self.quat_arr[i, :])
            omega_mb = self.omega_mb_arr[i, :]
            # logger.info("Derotating...")
            logger.debug("quat:     {}".format(quat))
            logger.debug("omega_mb: {}".format(omega_mb))
            logger.debug("pivot:    {}".format(self.pivot))
            snap['pos'], snap['vel'] = derotate_pos_and_vel(snap['pos'], snap['vel'], quat, omega_mb, self.pivot)

        if on_orbit_plane:
            logger.debug("Rotating on the plane of the orbit...")
            snap['pos'], snap['vel'] = rotate_on_orbit_plane(snap['pos'], snap['vel'])



def specific_angmom(snap):
    return pynbody.array.SimArray(pynbody.analysis.angmom.ang_mom_vec(snap), units=gadget_angmom_units)/snap['mass'].sum()


def compute_angmom(sim, derotator=None, on_orbit_plane=True, radius=10, initial_rotation_simname=""):
    """
    Returns the angular momentum of stars and gas inside a sphere of `radius` center in the center of the stars.
    Units are: 10**10*u.solMass * u.km/u.s * u.kpc
    Side effect: the simulation snap_list will be full of None
    """

    sphere = pynbody.filt.Sphere(radius * pynbody.units.kpc)
    angmom = defaultdict(list)

    families =  ('s', 'g', 'dm')

    # if derotator is not None:
    #     logger.info("Initialize derotator...")
    #     derotator = Derotator(sim)
    #     # derotator.derotate_sim(sim, on_orbit_plane=on_orbit_plane)

    logger.info('Computing angmom...')
    # From this: https://stackoverflow.com/a/42731787
    # sl = sim._snap_indexes
    # snap_indexes = list(range(sl.start or 0, sl.stop or len(sim), sl.step or 1))
    for i, snap in enumerate(tqdm.tqdm(sim)):
        try:
            if derotator is not None:
                snap['pos'], snap['vel'] = derotator.derotate_snap(snap, i)

            # Here I need to center on velocity
            pynbody.analysis.halo.center(snap.s)

            angmom['Ltot' + '_c'].append(pynbody.analysis.angmom.ang_mom_vec(snap[sphere]))
            angmom['j' + '_c'].append(specific_angmom(snap[sphere]))
            angmom['j'].append(specific_angmom(snap))

            for f in families:
                fam = pynbody.family.get_family(f)
                angmom['L'+ f + '_c'].append(pynbody.analysis.angmom.ang_mom_vec(snap[fam][sphere]))
                angmom['j'+ f + '_c'].append(specific_angmom(snap[fam][sphere]))
                angmom['j'+ f].append(specific_angmom(snap[fam]))

            baryon_snap = snap.g.union(snap.s)
            angmom['jb'].append(specific_angmom(baryon_snap))
            baryon_snap_sphere = baryon_snap[sphere]
            angmom['Lb' + '_c'].append(pynbody.analysis.angmom.ang_mom_vec(baryon_snap_sphere))
            angmom['jb' + '_c'].append(specific_angmom(baryon_snap_sphere))

        except ValueError as e:
            # Usually is 'Insufficient particles around center to get velocity'
            logger.warning(f"{i} {e}")
            # Get at least the time
            for k in angmom.keys():
                angmom[k].append([np.nan] * 3)
        del snap
        sim.snap_list[i] = None  # destroying references to the snap and the list
        if i % 10 == 0:
            gc.collect()

    for k, v in angmom.items():
        angmom[k] = np.vstack(v)

    # L, Lg = np.vstack(angmom), np.vstack(angmom_g)
    # print("before derotation")
    # print(L, Lg)
    # This is important in order to compare angular momentum from Moria to MovingBox
    # It is important to do this before the on_orbit_plane rotation
    if initial_rotation_simname:
        logger.info('Apply initial rotation...')
        quat_vp0 = get_initial_rotation(initial_rotation_simname)
        print(quat_vp0)
        for k, v in angmom.items():
            angmom[k] = rotate_vec(v, quat_vp0)
        # L = rotate_vec(L, quat_vp0)
        # Lg = rotate_vec(Lg, quat_vp0)

    if on_orbit_plane:
        logger.info('Rotating on orbit plane...')
        for k, v in angmom.items():
            angmom[k] = rotate_vec_on_orbit_plane(v)
        # L = rotate_vec_on_orbit_plane(L)
        # Lg = rotate_vec_on_orbit_plane(Lg)
    # print("After derotation")
    # print(L, Lg)

    return angmom
