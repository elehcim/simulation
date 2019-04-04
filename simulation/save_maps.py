import argparse
import functools
import os
import sys
import tqdm
import gc
from pprint import pprint

import matplotlib.pyplot as plt
import numpy as np
import pynbody
from astropy.table import Table, Column
from astropy import units as u

from simulation.luminosity import surface_brightness, kpc2pix, pix2kpc
from simulation.util import setup_logger, get_sim_name, to_astropy_quantity
from simulation.angmom import faceon, sideon

R_EFF_BORDER = 10

logger = setup_logger('__name__', logger_level='INFO')

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




class Snap:
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
            logger.info("Derotating...")
            logger.info("omega: {}".format(derot_param['omega']))
            logger.info("pivot: {}".format(derot_param['pivot']))
            s['vel'] -= np.cross(derot_param['omega'], s['pos'] - derot_param['pivot'])

        logger.info("Centering on stars")
        # vcen = pynbody.analysis.halo.vel_center(s, retcen=True)
        # logger.info("Original velocity center:", vcen)

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


def get_outname(snap_name, out_dir, band, width, resolution, suffix=None, stem_out = 'maps_'):
    out_name = stem_out + os.path.basename(snap_name) + '_{}_w{}_r{}'.format(band, width, resolution)
    out_name = os.path.join(os.path.expanduser(out_dir), out_name)
    if not os.path.isdir(out_name):
        logger.info('Creating folder {}'.format(out_name))
        os.makedirs(out_name, exist_ok=True)

    if suffix:
        out_name += suffix
    return out_name


def single_snap_maps(snap_name, width, resolution, band='v', side=True, face=None, omega=None, pivot=None,
                     sb_range=(18, 29), v_los_range=(-15, 15), sigma_range=(10, 40)):

    if omega is None or pivot is None:
        derot_param = None
    else:
        derot_param = {'omega': omega, 'pivot': pivot}

    snap = Snap(os.path.expanduser(snap_name), sphere_edge=R_EFF_BORDER, derot_param=derot_param)

    if side and face:
        print("Option 'side' and 'face' are mutually exclusive", file=sys.stderr)
        sys.exit(2)
    elif side:
        snap.sideon()
    elif face:
        snap.faceon()

    im = Imaging(snap.subsnap, width=width, resolution=resolution)

    return im


COLUMNS_UNITS = dict(vlos=u.km/u.s, sig=u.km/u.s, mag=u.mag * u.arcsec**-2, lum=u.solLum * u.pc**-2)


def simulation_maps(sim_path, width, resolution, band='v', side=True, face=None, omega_dir=None, pivot=None,
                     sb_range=(18, 29), v_los_range=(-15, 15), sigma_range=(10, 40)):

    if not os.path.isdir(sim_path):
        raise RuntimeError('Simulation path should be a directory')

    from simulation.util import snapshot_list
    snap_list = snapshot_list(sim_path, include_dir=True)

    if os.path.isdir(os.path.expanduser(omega_dir)):
        omega_dir = os.path.expanduser(omega_dir)
        omega_file = os.path.join(omega_dir, get_sim_name(sim_path)+'_omega.fits')
    else:
        omega_file = None

    if os.path.isfile(omega_file):
        logger.info('Reading omega table: {}'.format(omega_file))
        omega_arr = Table.read(omega_file)['omega'].data
        assert len(omega_arr) == len(snap_list)
    else:
        logger.warning('Cannot find omega table, not derotating...')
        omega_arr = None
    if pivot is not None:
        pivot = np.array(pivot.split(), dtype=np.float64)
    else:
        logger.info('No pivot provided, not derotating...')
        pivot = None

    nan_arr = np.empty((resolution, resolution), dtype=np.float32)
    nan_arr[:] = np.nan

    sim_name = get_sim_name(sim_path)

    maps_dict = dict(vlos=list(), sig=list(), mag=list(), lum=list())

    data_out_name = get_outname('data', out_dir=sim_name, band=band, width=width, resolution=resolution)
    print(data_out_name)

    for i, snap_name in enumerate(tqdm.tqdm(snap_list[:4])):

        snap = pynbody.load(snap_name)
        time = snap.header.time
        if omega_arr is not None:
            omega = omega_arr[i, :]
        else:
            omega = None

        try:

            im = single_snap_maps(snap_name=snap_name,
                                  width=width,
                                  resolution=resolution,
                                  band=band,
                                  side=side,
                                  face=face,
                                  sb_range=sb_range,
                                  v_los_range=v_los_range,
                                  sigma_range=sigma_range,
                                  omega=omega,
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
    group.add_argument("--snap", dest='snap_name', help="Path to the simulation snapshot")
    group.add_argument("--sim", "-s", dest='sim_path', help="Path to the simulation snapshot")
    parser.add_argument("--width", '-w', default=10, type=float, help='In kpc')
    parser.add_argument("--resolution", '-r', default=400, type=int)
    parser.add_argument("--band", "-b", default='v')
    parser.add_argument("--omega-dir", default='~/sim/analysis/ng_ana/data/omega', help='Directory of precomputed moving boxes omegas')
    parser.add_argument('--omega', help='Omega value in code units (space separated, e.g. "0 0 1.2")', type=str, default=None)
    parser.add_argument('--pivot', help='Coordinates of the pivot point (space separated, e.g. "30 30 30")', type=str, default=None)
    parser.add_argument("--out-dir", default=None)
    parser.add_argument('--side', action='store_true')
    parser.add_argument('--face', action='store_true')
    parser.add_argument("--sb-range", default=(18, 29), type=tuple)
    parser.add_argument("--v-los-range", default=(-15, 15), type=tuple)
    parser.add_argument("--sigma-range", default=(10, 40), type=tuple)


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
                              omega=np.array(args.omega.split(), dtype=np.float64),
                              pivot=np.array(args.pivot.split(), dtype=np.float64),
                              sb_range=args.sb_range,
                              v_los_range=args.v_los_range,
                              sigma_range=args.sigma_range,
                         )
    elif args.sim_path is not None:
        simulation_maps(sim_path=args.sim_path,
                        width=args.width,
                        resolution=args.resolution,
                        band=args.band,
                        side=args.side,
                        face=args.face,
                        omega_dir=args.omega_dir,
                        pivot=args.pivot,
                        sb_range=args.sb_range,
                        v_los_range=args.v_los_range,
                        sigma_range=args.sigma_range,
                        )

    return im

if __name__ == '__main__':
    main()