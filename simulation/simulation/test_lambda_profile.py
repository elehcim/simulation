from lambda_r import ss_angmom_profile, ss_angmom, integrate_annulus
from lambda_r import create_apertures
from photutils import aperture_photometry

import matplotlib.pyplot as plt
from astropy.table import Table
import os
import astropy.units as u
import numpy as np
from util import to_astropy_quantity

resolution = 200
n_annuli=20
a_min=10
ellip = 0.2
theta=0
center=(0,0)
bit_less_than_half_diagonal = int(resolution*1.3)

smajax = np.linspace(a_min, bit_less_than_half_diagonal, n_annuli)


tbl = Table.read('/home/michele/sim/analysis/second_ssam_run/m69p3/maps_data_v_w10_r200_a20.fits_img.fits.gz')

v_los_map = tbl['vlos'].data[1]
v_disp_map= tbl['sig'].data[1]
sb_mag= tbl['mag'].data[1]


lum = u.Quantity(sb_mag, unit='mag/arcsec**2')
v_los = u.Quantity(v_los_map, unit = 'km/s')
v_disp = u.Quantity(v_disp_map, unit = 'km/s')

apertures = create_apertures(center, smajax, ellip, theta)
flux_table = aperture_photometry(lum, apertures)



lum_annuli = integrate_annulus(lum, center, smajax, ellip, theta)
v_los_annuli = integrate_annulus(v_los, center, smajax, ellip, theta)
v_disp_annuli = integrate_annulus(v_disp, center, smajax, ellip, theta)

l_prof = ss_angmom_profile(lum_annuli, smajax, v_los_annuli, v_disp_annuli)
l = ss_angmom(lum_annuli, smajax, v_los_annuli, v_disp_annuli)

plt.plot(smajax, l_prof)