import argparse
import functools
import os
import sys
from copy import deepcopy
from pprint import pprint
import tqdm

import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, Column
from photutils import EllipticalAperture, data_properties

from simulation.lambda_r_photometric import print_fit_results, fit_sersic_2D, create_apertures
from simulation.util import setup_logger

logger = setup_logger('__name__', logger_level='INFO')


class Photometry:
    sersic1D = None
    sersic2D = None

    def __init__(self, band, resolution, n_annuli=20, a_min=10, N_0=1, ELLIP_0=0, THETA_0=0):
        """
        a_min
        in pixels the minimum of the semimajor axis and the interval between one aperture and the other.
        For now it is linked with the resolution.
        """
        self.band = band
        self.n_annuli = n_annuli
        self.resolution = resolution
        self._n_0 = N_0
        self._ellip_0 = ELLIP_0
        self._theta_0 = THETA_0
        bit_less_than_half_diagonal = int(resolution*1.3)
        self.smajax = np.linspace(a_min, bit_less_than_half_diagonal, n_annuli)
        self._sersic2D = None


    @property
    def data_properties(self):
        return
    def fit(self, sb_lum, r_eff_pix, fit_profile=False, snap=None, fixed=None):
        if fit_profile:
            sersic1D = self.fit_profile(snap, r_eff_pix)
            self._n_0 = sersic1D.n
        logger.info("Fitting Sersic2D")
        sersic2D = fit_sersic_2D(sb_lum, r_eff=r_eff_pix, n=self._n_0, resolution=self.resolution,
                                 ellip=self._ellip_0, theta=self._theta_0, fixed=fixed)

        self._sersic2D = sersic2D

        print_fit_results(sersic2D)

        if 0 <= sersic2D.ellip.value <= 1:
            self.ellip = sersic2D.ellip.value
            self.theta = sersic2D.theta.value
        elif sersic2D.ellip.value < 0:
            logger.warning("ellipticity < 0: swapping minor <-> major axis")
            self.ellip = 1 / (1 - sersic2D.ellip.value)
            self.theta = np.pi / 2 + sersic2D.theta.value
        else:
            logger.warning("ellipticity > 1: swapping minor <-> major axis")
            self.ellip = 1 / sersic2D.ellip.value
            self.theta = np.pi / 2 + sersic2D.theta.value
        # Use half of the quadrant.
        self.theta = self.theta % np.pi
        self.center = (sersic2D.x_0.value, sersic2D.y_0.value)
        self.n = sersic2D.n.value
        return sersic2D

    # FIXME
    def fit_profile(self, snap, r_eff_pix):
        raise NotImplementedError
        # logger.info("Creating profile")
        # r_bins, sbp = sb_profile(snap, band=self.band)
        # if SHOW:
        #     plt.plot(r_bins, sbp, linewidth=2);
        #     plt.ylabel("I [$L_{{\odot,{}}}/pc^2$]".format(self.band))
        #     plt.xlabel("r/kpc")
        #     plt.title('Surface brightness profile')
        # logger.info("Fitting Sersic1D")
        # # TODO use r_eff_3D here?
        # sersic1D = fit_sersic_1D(r_eff_kpc, r_bins, sbp, show=SHOW)

        # self.sersic1D = sersic1D

        # print_fit_results(sersic1D)
        # return sersic1D

    @property
    def apertures(self):
        apertures = create_apertures(self.center, self.smajax, self.ellip, self.theta)
        return apertures

    def get_params(self):
        return self.center, self.smajax, self.ellip, self.theta, self.n



def parse_args(cli=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--maps", help="Path to the maps")


    args = parser.parse_args(cli)
    return args


def main(cli=None):
    args = parse_args(cli)

    pprint(args.__dict__)

    tbl = Table.read(args.maps)
    prop_list = list()
    for i, img in enumerate(tqdm.tqdm(tbl['lum'])):
        cat_tbl = data_properties(img).to_table()
        cat_tbl['id'] = i + 1
        prop_list.append(cat_tbl)
    prop = vstack(prop_list)
    return prop

if __name__ == '__main__':
    prop = main()