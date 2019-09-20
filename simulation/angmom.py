import pynbody
from pynbody.analysis.angmom import calc_faceon_matrix, calc_sideon_matrix

from astropy import units as u
from simulation.util import get_quat_omega_pivot, get_sim_name, setup_logger
from simulation.derotate_simulation import (derotate_pos_and_vel, derotate_simulation,
                                           rotate_vec_on_orbit_plane, rotate_on_orbit_plane)
from astropy.table import Table
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
    def __init__(self, sim_name, slicer=slice(None)):
        quat_arr, omega_mb_arr, pivot = get_quat_omega_pivot(sim_name)
        self.quat_arr = quat_arr[slicer]
        self.omega_mb_arr = omega_mb_arr[slicer]
        self.pivot = pivot


    def derotate_snap(self, snap):
        snap['pos'], snap['vel'] = derotate_pos_and_vel(snap['pos'], snap['vel'], self.quat, self.omega_mb, self.pivot)
        return snap

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
            logger.info("Rotating on the plane of the orbit...")
            snap['pos'], snap['vel'] = rotate_on_orbit_plane(snap['pos'], snap['vel'])



def compute_angmom(sim, slicer=slice(None), derotate=True, on_orbit_plane=True, radius=10):
    """
    Returns the angular momentum of stars and gas inside a sphere of `radius` center in the center of the stars.
    Units are: 10**10*u.solMass * u.km/u.s * u.kpc
    Side effect: the simulation snap_list will be full of None
    """

    sim_name = get_sim_name(sim.sim_id)
    sphere = pynbody.filt.Sphere(radius * pynbody.units.kpc)
    angmom = list()
    angmom_g = list()

    if derotate:
        derotator = Derotator(sim_name, slicer=slicer)
        derotator.derotate_sim(sim, on_orbit_plane=on_orbit_plane)

    logger.info('Computing angmom...')
    for i, snap in enumerate(tqdm.tqdm(sim)):
        try:
            # Here I need to center on velocity
            pynbody.analysis.halo.center(snap.s)
            # sideon(snap.s[sphere])
            angmom.append(pynbody.analysis.angmom.ang_mom_vec(snap.s[sphere]))

            # sideon(snap.s[sphere])
            angmom_g.append(pynbody.analysis.angmom.ang_mom_vec(snap.g[sphere]))

        except ValueError as e:
            # Usually is 'Insufficient particles around center to get velocity'
            print(e)
            # Get at least the time
            angmom.append([np.nan] * 3)
            angmom_g.append([np.nan] * 3)
        del snap
        sim.snap_list[i] = None  # destroying references to the snap and the list
        if i % 10 == 0:
            gc.collect()

    L, Lg = np.vstack(angmom), np.vstack(angmom_g)

    if on_orbit_plane:
        L = rotate_vec_on_orbit_plane(L)
        Lg = rotate_vec_on_orbit_plane(Lg)

    return L, Lg