import argparse
import functools
import os
import sys
import tqdm
import gc
from pprint import pprint

import matplotlib.pyplot as plt
import numpy as np
import quaternion
import pynbody
from astropy.table import Table, Column
from astropy import units as u

from simulation.luminosity import surface_brightness, kpc2pix, pix2kpc
from simulation.util import setup_logger, get_sim_name, to_astropy_quantity, get_pivot, get_quat
from simulation.angmom import faceon, sideon

R_EFF_BORDER = 10

logger = setup_logger('save_maps', logger_level='INFO')

class Imaging:
    def __init__(self, snap, width, resolution):
        self.width = width
        self.resolution = resolution
        self._snap = snap

    @functools.lru_cache(1)
    def v_los_map(self):
        return pynbody.plot.sph.image(self._snap.s, qty='vz', av_z=True, width=self.width,
                                      resolution=self.resolution, noplot=True, log=False)

    @functools.lru_cache(1)
    def v_disp_map(self):
        return pynbody.plot.sph.image(self._snap.s, qty='v_disp', av_z=True, width=self.width,
                                      resolution=self.resolution, noplot=True, log=False)

    @functools.lru_cache(1)
    def sb_lum(self):
        return surface_brightness(self._snap.s, width=self.width, resolution=self.resolution,
                                  lum_pc2=True, noplot=True)

    @functools.lru_cache(1)
    def sb_mag(self):
        return surface_brightness(self._snap.s, width=self.width, resolution=self.resolution,
                                  lum_pc2=False, noplot=True)

    @classmethod
    def from_fits(cls):
        pass

def rotate_vec(vec, quat):
    """Rotate a numpy array of 3-vectors `vec` given a quaternion `quat`"""
    v = np.hstack([np.zeros((len(vec), 1)), vec])
    vq = quaternion.as_quat_array(v)
    new_vec = quat * vq * quat.conj()

    # remove w component
    return quaternion.as_float_array(new_vec)[:, 1:]


def derotate_snap():
    pass

class Snap:
    """
    The main reason for this class is to have a uniform preprocessing.
    In particular the selection and the orientation should be well defined
    """
    def __init__(self, snap_name, sphere_edge, derot_param=None):
        logger.info("Opening file {}".format(snap_name))
        s = pynbody.load(snap_name)
        self._snap = s
        # max_boxsize = 4000
        # s.properties['boxsize'] = pynbody.units.Unit("{} kpc".format(max_boxsize))
        # s.physical_units()
        self.time = s.header.time
        self.time_gyr = s.properties['time'].in_units('Gyr')
        logger.info("{:.2f} Gyr".format(self.time))

        if derot_param is not None:
            quat = np.quaternion(*derot_param['quat'])
            logger.info("Derotating...")

            logger.info("quat: {}".format(quat))
            logger.info("pivot: {}".format(derot_param['pivot']))

            s['pos'] = rotate_vec(s['pos'] - derot_param['pivot'], quat)
            s['vel'] = rotate_vec(s['vel'], quat)

            # s['vel'] -= np.cross(derot_param['omega'], s['pos'] - derot_param['pivot'])

        logger.info("Centering on stars")
        # vcen = pynbody.analysis.halo.vel_center(s, retcen=True)
        # logger.info("Original velocity center:", vcen)

        # FIXME Maybe this is not necessary
        pynbody.analysis.halo.center(s.s, vel=True)

        vcen_new = pynbody.analysis.halo.vel_center(s.s, retcen=True)
        logger.info("New velocity center: {}".format(vcen_new))

        # self.subsnap = s[pynbody.filt.Cuboid('{} kpc'.format(-cuboid_edge))]
        self.subsnap = s[pynbody.filt.Sphere('{} kpc'.format(sphere_edge))]


    def sideon(self):
        logger.info("Rotating sideon")
        # pynbody.analysis.angmom.sideon(self.subsnap.s, disk_size=self.subsnap.s['pos'].max())
        sideon(self.subsnap.s)

    def faceon(self):
        logger.info("Rotating faceon")
        # pynbody.analysis.angmom.faceon(self.subsnap.s, disk_size=self.subsnap.s['pos'].max())
        faceon(self.subsnap.s)

    @property
    def angmom(self):
        logger.info("Computing overall angular momentum")
        return pynbody.analysis.angmom.ang_mom_vec(self.subsnap)

    @property
    @functools.lru_cache(1)
    def r_eff_kpc3d(self):
        return pynbody.analysis.luminosity.half_light_r(self.subsnap.s, cylindrical=False)

    @property
    @functools.lru_cache(1)
    def r_eff_kpc(self):
        return pynbody.analysis.luminosity.half_light_r(self.subsnap.s, cylindrical=True)

    def magnitude(self, band):
        return pynbody.analysis.luminosity.halo_mag(self.subsnap.s, band=band)


def get_outname(snap_name, out_dir, band, width, resolution, suffix=None, stem_out='maps_'):
    out_name = stem_out + os.path.basename(snap_name) + '_{}_w{}_r{}'.format(band, width, resolution)
    out_name = os.path.join(os.path.expanduser(out_dir), out_name)
    if not os.path.isdir(out_name):
        logger.info('Creating folder {}'.format(out_name))
        os.makedirs(out_name, exist_ok=True)

    if suffix:
        out_name += suffix
    return out_name


def single_snap_maps(snap_name, width, resolution, band='v', side=True, face=None, rotation=None, quat=None, pivot=None):

    """
    Parameters
    ----------
    rotation: pynbody.transformation
        The rotation
    quat: np.array
    pivot: np.array

    """
    if quat is None or pivot is None:
        derot_param = None
    else:
        derot_param = {'quat': quat, 'pivot': pivot}

    snap = Snap(os.path.expanduser(snap_name), sphere_edge=R_EFF_BORDER, derot_param=derot_param)

    if rotation is not None:
        rotation.apply(snap)
    else:
        if side:
            snap.sideon()
        elif face:
            snap.faceon()
        else:
            logger.warning('No sideon or faceon indications: not rotating snap')

    im = Imaging(snap.subsnap, width=width, resolution=resolution)

    return im


COLUMNS_UNITS = dict(vlos=u.km/u.s, sig=u.km/u.s, mag=u.mag * u.arcsec**-2, lum=u.solLum * u.pc**-2)


def simulation_maps(sim_path, width, resolution, band='v', side=True, face=None, quat_dir=None, pivot=None):

    if not os.path.isdir(sim_path):
        raise RuntimeError('Simulation path should be a directory')

    from simulation.util import snapshot_list
    snap_list = snapshot_list(sim_path, include_dir=True)

    sim_name = get_sim_name(sim_path)

    quat_arr = get_quat(sim_name, quat_dir)
    if quat_arr is None:
        raise RuntimeError('Cannot find quaternion. First save it using save_quat.py')

    assert len(quat_arr) == len(snap_list)

    if pivot is None:
        pivot = get_pivot(sim_name)

    nan_arr = np.empty((resolution, resolution), dtype=np.float32)
    nan_arr[:] = np.nan


    maps_dict = dict(vlos=list(), sig=list(), mag=list(), lum=list())

    data_out_name = get_outname('data', out_dir=sim_name, band=band, width=width, resolution=resolution)
    print(data_out_name)

    for i, snap_name in enumerate(tqdm.tqdm(snap_list)):

        snap = pynbody.load(snap_name)
        time = snap.header.time
        if quat_arr is not None:
            quat = quat_arr[i, :]
        else:
            quat = None

        try:

            im = single_snap_maps(snap_name=snap_name,
                                  width=width,
                                  resolution=resolution,
                                  band=band,
                                  side=side,
                                  face=face,
                                  quat=quat,
                                  pivot=pivot,
                                  )
            # vlos = to_astropy_quantity(im.v_los_map())
            # sig = to_astropy_quantity(im.v_disp_map())
            # mag = im.sb_mag() * u.mag * u.arcsec**-2
            # lum = im.sb_lum() * u.solLum * u.pc**-2
            vlos = im.v_los_map()
            sig = im.v_disp_map()
            mag = im.sb_mag()
            lum = im.sb_lum()
            # Map saving
            maps_dict['vlos'].append(vlos)
            maps_dict['sig'].append(sig)
            maps_dict['mag'].append(mag)
            maps_dict['lum'].append(lum)
            mtbl_loc = Table([vlos, sig, mag, lum], names=['vlos','sig','mag', 'lum'], meta={'time':time, 'time_u': str(u.kpc*u.s*u.km**-1)})
            for col_name, col in mtbl_loc.columns.items():
                col.unit = COLUMNS_UNITS[col_name]

            mtbl_loc.write(os.path.join(data_out_name, 'maps{}_img.fits.gz'.format(snap_name[-4:])), overwrite=True)

        except ValueError as e:
            # Usually is 'Insufficient particles around center to get velocity'
            logger.error(snap_name)
            logger.error(e)

            for k, v in maps_dict.items():
                v.append(nan_arr)

        del snap
        if i % 5 == 0:
            gc.collect()

    fits_data_file = data_out_name+'_img.fits'
    logger.info("Writing fits with all the maps... ({}) ".format(fits_data_file))
    mtbl = Table(maps_dict)
    for col_name, col in mtbl.columns.items():
        col.unit = COLUMNS_UNITS[col_name]
    mtbl.write(fits_data_file, overwrite=True)


def parse_args(cli=None):
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    angmom_group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument("--snap", dest='snap_name', help="Path to the single snapshot to analyze")
    group.add_argument("--sim", "-s", dest='sim_path', help="Path to the simulation snapshot")
    parser.add_argument("--width", '-w', default=10, type=float, help='In kpc')
    parser.add_argument("--resolution", '-r', default=400, type=int)
    parser.add_argument("--band", "-b", default='v')
    parser.add_argument("--quat-dir", default='~/sim/analysis/ng_ana/data/quat', help='Directory of precomputed moving boxes quaternions')
    parser.add_argument('--quat', help='Quaternion value (space separated, e.g. "1 0 0 1.2")', type=str, default=None)
    parser.add_argument('--pivot', help='Coordinates of the pivot point (space separated, e.g. "30 30 30")', type=str, default=None)
    parser.add_argument("--out-dir", default=None)
    angmom_group.add_argument('--side', action='store_true')
    angmom_group.add_argument('--face', action='store_true')

    args = parser.parse_args(cli)
    return args


def main(cli=None):
    args = parse_args(cli)

    pprint(args.__dict__)
    im = None

    if args.snap_name is not None:
        im = single_snap_maps(snap_name=args.snap_name,
                              width=args.width,
                              resolution=args.resolution,
                              band=args.band,
                              side=args.side,
                              face=args.face,
                              quat=np.array(args.quat.split(), dtype=np.float64),
                              pivot=np.array(args.pivot.split(), dtype=np.float64),
                              )
    elif args.sim_path is not None:
        simulation_maps(sim_path=args.sim_path,
                        width=args.width,
                        resolution=args.resolution,
                        band=args.band,
                        side=args.side,
                        face=args.face,
                        quat_dir=args.quat_dir,
                        pivot=args.pivot,
                        )

    return im

if __name__ == '__main__':
    main()