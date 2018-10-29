import argparse
import functools
import logging
import os
import sys
from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np
import pynbody
from photutils import EllipticalAnnulus

from lambda_r import print_fit_results, plot_angmom, compute_stellar_specific_angmom, plot_maps, fit_sersic_2D, \
    fit_sersic_1D, sb_profile
from luminosity import surface_brightness, kpc2pix, pix2kpc

logFormatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s]  %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logger.addHandler(consoleHandler)


class Photometry:
    sersic1D = None
    sersic2D = None

    def __init__(self, band, a_delta, a_min=30, a_max=200, N_0=1, ELLIP_0=0, THETA_0=0):
        self.band = band
        self.a_delta = a_delta
        self._n_0 = N_0
        self._ellip_0 = ELLIP_0
        self._theta_0 = THETA_0
        self.smajax = np.arange(a_min, a_max, a_delta)

    def fit(self, sb_lum, r_eff_pix, resolution, fit_profile=False, snap=None):
        if fit_profile:
            sersic1D = self.fit_profile(snap, r_eff_pix)
            self._n_0 = sersic1D.n
        logger.info("Fitting Sersic2D")
        sersic2D = fit_sersic_2D(sb_lum, r_eff=r_eff_pix, n=self._n_0, resolution=resolution, 
            ellip=self._ellip_0, theta=self._theta_0)

        self._sersic2D = sersic2D

        print_fit_results(sersic2D)

        if 0 <= sersic2D.ellip.value <= 1:
            self.ellip = sersic2D.ellip.value
            self.theta = sersic2D.theta.value
        elif sersic2D.ellip.value < 0:
            logger.warning("ellipticity < 0: swapping minor <-> major axis")
            self.ellip = 1/(1-sersic2D.ellip.value)
            self.theta = np.pi/2 + sersic2D.theta.value
        else:
            logger.warning("ellipticity > 1: swapping minor <-> major axis")
            self.ellip = 1/sersic2D.ellip.value
            self.theta = np.pi/2 + sersic2D.theta.value
        self.center = (sersic2D.x_0.value, sersic2D.y_0.value)

        return sersic2D

    # FIXME
    def fit_profile(self, snap, r_eff_pix):
        logger.info("Creating profile")
        r_bins, sbp = sb_profile(snap, band=self.band)
        if SHOW:
            plt.plot(r_bins, sbp, linewidth=2);
            plt.ylabel("I [$L_{{\odot,{}}}/pc^2$]".format(self.band))
            plt.xlabel("r/kpc")
            plt.title('Surface brightness profile')
        logger.info("Fitting Sersic1D")
        # TODO use r_eff_3D here?
        sersic1D = fit_sersic_1D(r_eff_kpc, r_bins, sbp, show=SHOW)

        self.sersic1D = sersic1D

        print_fit_results(sersic1D)
        return sersic1D

    @property
    def apertures(self):
        sminax = self.smajax * np.sqrt(1 - self.ellip)
        apertures = list()
        for a, b in zip(self.smajax, sminax):
            apertures.append(EllipticalAnnulus(self.center,
                                       a_in=a-self.a_delta, a_out=a+self.a_delta,
                                       b_out=b, theta=self.theta))
        return apertures

    def get_params(self):
        return self.center, self.smajax, self.ellip, self.theta


class Imaging:
    def __init__(self, snap, width, resolution):
        self.width = width
        self.resolution = resolution
        self._snap = snap

    def v_los_map(self):
        return pynbody.plot.sph.image(self._snap.s, qty='vz', av_z=True, width=self.width,
            resolution=self.resolution, noplot=True, log=False)

    def v_disp_map(self):
        return pynbody.plot.sph.image(self._snap.s, qty='v_disp', av_z=True, width=self.width,
            resolution=self.resolution, noplot=True, log=False)

    def sb_lum(self):
        return surface_brightness(self._snap.s, width=self.width, resolution=self.resolution,
            lum_pc2=True, noplot=True)

    def sb_mag(self):
        return surface_brightness(self._snap.s, width=self.width, resolution=self.resolution,
            lum_pc2=False, noplot=True)

class Snap:
    def __init__(self, snap_name, cuboid_edge):
        logger.info("Opening file {}".format(snap_name))
        s = pynbody.load(snap_name)
        self._snap = s
        max_boxsize = 4000
        s.properties['boxsize'] = pynbody.units.Unit("{} kpc".format(max_boxsize))
        s.physical_units()
        self.time = s.properties['time'].in_units('Gyr')
        logger.info("{:.2f} Gyr".format(self.time))

        pynbody.analysis.halo.center(s.s)#, vel=False)
        self.subsnap = s[pynbody.filt.Cuboid('{} kpc'.format(-cuboid_edge))] # -width*1.1

    def sideon(self):
        logger.info("Rotating sideon")
        pynbody.analysis.angmom.sideon(self.subsnap.s, disk_size=self.subsnap.s['pos'].max())

    def faceon(self):
        logger.info("Rotating faceon")
        pynbody.analysis.angmom.faceon(self.subsnap.s, disk_size=self.subsnap.s['pos'].max())

    @property    
    def angmom(self):
        logger.info("Computing overall angular momentum")
        return pynbody.analysis.angmom.ang_mom_vec(self.subsnap)

    @property
    @functools.lru_cache(1)
    def r_eff_kpc3d(self):
        return pynbody.analysis.luminosity.half_light_r(self.subsnap, cylindrical=False)

    @property
    @functools.lru_cache(1)
    def r_eff_kpc(self):
        return pynbody.analysis.luminosity.half_light_r(self.subsnap, cylindrical=True)


class SSAM:
    def __init__(self, snap, photometry, imaging):
        self.photometry = photometry
        self.imaging = imaging
        self.w = self.imaging.width
        self.res = self.imaging.resolution
        self.snap = snap
        self.sb_mag = self.imaging.sb_mag()
        self.v_los_map, self.v_disp_map = self.imaging.v_los_map(), self.imaging.v_disp_map()
        self.time = snap.time
        self.lambda_R = None

    def compute_lambda(self, fit_profile=False):
        sb_lum = self.imaging.sb_lum()
        logger.info("Computed R_eff:")
        logger.info(" 2D: {:.4f} kpc".format(self.snap.r_eff_kpc))
        # logger.info(" 3D: {:.4f} kpc".format(self.snap.r_eff_kpc3d))
        
        r_eff_pix = kpc2pix(self.snap.r_eff_kpc,  # use partial?
                            width=self.w,
                            resolution=self.res)
        self.photometry.fit(sb_lum,
                            r_eff_pix=r_eff_pix,
                            resolution=self.res,
                            fit_profile=fit_profile,
                            snap=self.snap.subsnap)
        
        center, smajax, ellip, theta = self.photometry.get_params()
        lambda_R = compute_stellar_specific_angmom(center, self.sb_mag, self.v_los_map, self.v_disp_map, 
                                                   smajax, ellip, self.photometry.a_delta, theta)
        self.lambda_R = lambda_R
        logger.info("Lambda_R = {:.4f}".format(lambda_R))
        return lambda_R

    def plot_maps(self, save_fig=None, sb_range=None, v_los_range=None, sigma_range=None):
        grid = plot_maps(self.sb_mag, self.v_los_map, self.v_disp_map, 
                         self.w, self.res,
                         sb_range, v_los_range, sigma_range)
        annuli_to_plot = deepcopy(self.photometry.apertures)
        for aper in annuli_to_plot:
            for attr in ('a_in', 'a_out', "b_in", "b_out"):
                qty = getattr(aper, attr)
                setattr(aper, attr, pix2kpc(qty, width=self.w, resolution=self.res))
            aper.positions = np.array([[0,0]])
            aper.plot(color='white', ax=grid[0], alpha=0.2)

        plot_angmom(self.snap.subsnap.s, grid[0])
        fig = plt.gcf()
        fig.suptitle('t={:.2f} Gyr  $\lambda_R$={:.4f}'.format(self.time, self.lambda_R))
        plt.subplots_adjust(top=0.80)
        if save_fig:
            logger.info("Saving figure: {}".format(save_fig))
            plt.savefig(save_fig)
        else:
            plt.show()

if __name__ == '__main__':
    # snap_name = "/home/michele/sim/MoRIA/M1-10_Verbeke2017/M10sim41001/snapshot_0036"
    # snap_name = "/mnt/data/MoRIA/M1-10_Verbeke2017/M10sim41001/snapshot_0036"
    # snap_name = "/home/michele/sim/MySimulations/hi_osc/mb.69002_p200_a800_r600/out/snapshot_0039"
    # width=10
    # resolution = 500
    # a_delta = 20
    # band = 'v'

    parser = argparse.ArgumentParser()
    parser.add_argument("--snap", "-s", help="Path to the simulation snapshot", required=True)
    parser.add_argument("--width", '-w', default=10, type=float)
    parser.add_argument("--resolution", '-r', default=400, type=int)
    parser.add_argument("--band", "-b", default='v')
    parser.add_argument("--a-delta", default=20, type=int)
    parser.add_argument("--cuboid-factor", default=1.5, type=float)
    parser.add_argument("--out-dir", default=None)
    parser.add_argument('--side', action='store_true')
    parser.add_argument('--face', action='store_true')

    args = parser.parse_args()

    snap_name = args.snap
    width = args.width
    resolution = args.resolution
    a_delta = args.a_delta
    band = args.band
    out_dir = args.out_dir

    stem_out = 'maps_'
    out_name = stem_out + os.path.basename(snap_name) + '_{}_w{}_r{}_a{}'.format(band, width, resolution, a_delta)
    if out_dir is not None:
        if os.path.isdir(out_dir):
           out_name = os.path.join(os.path.expanduser(out_dir), out_name)

    snap = Snap(snap_name, cuboid_edge=width*1.1)

    if args.side and args.face:
        print("Option 'side' and 'face' are mutually exclusive", file=sys.stderr)
        sys.exit(2)
    elif args.side:
        snap.sideon()
    elif args.face:
        snap.faceon()

    im = Imaging(snap.subsnap, width=width, resolution=resolution)
    ph = Photometry(band=band, a_delta=a_delta)

    ssam = SSAM(snap, photometry=ph, imaging=im)

    ssam.compute_lambda(fit_profile=False)
    
    print(ssam.time, ssam.lambda_R)
    
    ssam.plot_maps(save_fig=out_name, sb_range=(18,29),
                                      v_los_range=(-15,15),
                                      sigma_range=(10,40))

