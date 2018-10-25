import pynbody
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.modeling import models, fitting
from photutils import isophote, aperture_photometry
from photutils import CircularAperture, EllipticalAperture, EllipticalAnnulus
from photutils.isophote import EllipseGeometry, Ellipse
from luminosity import surface_brightness, color_plot, kpc2pix

import logging
logFormatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s]  %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logger.addHandler(consoleHandler)


NP_ERRSTATE = {"divide":'ignore', "over":'ignore', "under":'ignore'}

N_0 = 1
ELLIP_0 = 0
THETA_0 = 0

FIT_PROFILE = True

SHOW = False

_sun_abs_mag = {'u':5.56,'b':5.45,'v':4.8,'r':4.46,'i':4.1,'j':3.66,
               'h':3.32,'k':3.28}
_arcsec2_over_pc2_at_10pc = (np.tan(np.pi/180/3600)*10.0)**2

def bounding_box(snap):
    unit = snap['pos'].units
    for coord in 'x y z'.split():
        print("{}: {:10.2f}, {:10.2f} ({})".format(coord, snap[coord].min(), snap[coord].max(), unit))
    return [(float(snap[coord].min()), float(snap[coord].max())) for coord in 'x y z'.split()]

def to_astropy_quantity(simarr, units=None):
    return u.Quantity(simarr.view(type=np.ndarray), unit=units if units is not None else str(simarr.units))

def plot_fit(img, fit):
    # Plot the data with the best-fit model
    res_x, res_y = img.shape
    y, x = np.mgrid[:res_x, :res_y]
    plt.figure(figsize=(12,2.5))
    plt.subplot(1,3,1)
    data = plt.imshow(img, origin='lower', interpolation='nearest')
    plt.title("Data")
    plt.colorbar(data)
    plt.subplot(1,3,2)
    mod = plt.imshow(fit(x, y), origin='lower', interpolation='nearest')
    plt.colorbar(mod)
    plt.title("Model")
    plt.subplot(1,3,3)
    residual = img - fit(x, y)
    res = plt.imshow(residual, origin='lower', interpolation='nearest')
    plt.colorbar(res)
    plt.title("Residual");
    print("residual: ({:.4f}, {:.4f})".format(residual.min(), residual.max()))

def plot_angmom(snap, ax):
    """Plot the projected angular momentum on the map on `ax`, after having
    normalized it to one. It means that short arrows plotted means L almost
    aligned with line-of-sight-direction."""
    L = pynbody.analysis.angmom.ang_mom_vec(snap)
    print("L:", L)
    norm = np.linalg.norm(L)
    ax.arrow(0, 0, L[0]/norm, L[1]/norm, head_width=0.2, head_length=.2, color='red');

def ss_angmom(flux, r, v_los, v_disp):
    return np.sum(flux * r * np.abs(v_los)) / np.sum(flux * r * np.sqrt(v_los**2 + v_disp**2))

# -2.5*np.log10(arcsec2_over_pc2_at_10pc) == 21.572

def sb_profile(band):
    binning='equaln' # contain equal numbers of particles
    sun_mag = _sun_abs_mag[band]
    ps = pynbody.analysis.profile.Profile(subsnap.s, type=binning, max=4*r_eff_kpc, bin=100)
    r = ps['rbins'].in_units('kpc')
    # sbp = 10**(0.4*(sun_abs_mag - 2.5*np.log10(arcsec2_over_pc2_at_10pc) - ps['sb,' + band] ))
    sbp = 10**(0.4*(sun_mag + 21.572 - ps['sb,' + band] ))
    plt.plot(r, sbp, linewidth=2);
    plt.ylabel("I [$L_{{\odot,{}}}/pc^2$]".format(band))
    plt.xlabel("r/kpc")
    plt.title('Surface brightness profile')
    return r, sbp


def fit_sersic_1D(r_eff_kpc3d, r, sbp, show=SHOW):
    s1D_init = models.Sersic1D(r_eff=r_eff_kpc3d, n=2, fixed={'amplitude':False, 'r_eff':False, 'n':False})
    # s1d_init = models.ExponentialCutoffPowerLaw1D(alpha=1, x_cutoff=3)
    fit_s1D = fitting.SLSQPLSQFitter()
    with np.errstate(**NP_ERRSTATE):
        sersic1D = fit_s1D(s1D_init, r, sbp)#, acc=1e-10)
    plt.plot(r, sbp, linewidth=2)
    plt.plot(r, sersic1D(r))
    if show:
        plt.show()
    return sersic1D


def fit_sersic_2D(sb, r_eff, n, resolution, ellip, theta, show=SHOW):
    y, x = np.mgrid[:resolution, :resolution]
    s_init = models.Sersic2D(r_eff=r_eff, n=n, x_0=resolution/2, y_0=resolution/2, ellip=ellip, theta=theta,
                             fixed={'amplitude':False, 'r_eff':True, 'n':False, 'x_0':True, 'y_0':True, 'ellip':False, 'theta':False})
                             # bounds={ 'ellip':(0,1)}) #, 'x_0':(200,300) 'y_0':(200,300),
    # s_init = models.Polynomial2D(degree=4)
    # s_init = models.Gaussian2D(x_mean=resolution/2, y_mean=resolution/2, theta=theta,
    #                            fixed={'amplitude':False, 'x_mean':False, 'y_mean':False, 'theta':False}, 
    #                            bounds={'theta':(0, np.pi), 'x_mean':(200,300), 'y_mean':(200,300)})

    fit_s = fitting.SLSQPLSQFitter()
    # notnans = np.isfinite(img)
    # sersic = fit_s(s_init, x[notnans], y[notnans], img[notnans])
    with np.errstate(**NP_ERRSTATE):
        sersic = fit_s(s_init, x, y, sb)#, weights=weights_norm, maxiter=200)
    print(sersic)
    if show:
        plot_fit(sb, sersic)
        plt.show()
    return sersic


def plot_aperture_geometry(sb, sersic, resolution, show=SHOW):
    if sersic.ellip.value <= 1:
        ellip = sersic.ellip.value
        theta = sersic.theta.value
    else:
        print("ellipticity > 1: swapping minor <-> major axis")
        ellip = 1/sersic.ellip.value
        theta = np.pi/2 + sersic.theta.value

    geometry = EllipseGeometry(x0=resolution/2, y0=resolution/2, sma=sersic.r_eff.value,
                               eps=ellip, pa=theta)
    # print(geometry)
    aper = EllipticalAperture((geometry.x0, geometry.y0), geometry.sma,
                               geometry.sma*(1 - geometry.eps),
                               geometry.pa)
    if show:
        plt.imshow(sb, origin='lower')
        aper.plot(color='white')
        plt.colorbar()
        plt.show()
    return aper

def integrate_annulus(qty, center, smajax, ellip, a_delta, theta):
    apertures = create_apertures(center, smajax, ellip, a_delta, theta)
    flux_table = aperture_photometry(qty, apertures)
#     for col in flux_table.colnames:
#         flux_table[col].info.format = '%.8g'  # for consistent table output
#     print(flux_table)
    return u.Quantity(np.array([flux_table['aperture_sum_{}'.format(i)][0].value for i in range(len(smajax))]), unit=flux_table['aperture_sum_0'].unit)

def create_apertures(center, smajax, ellip, a_delta, theta):
    if ellip > 1:
        print("ellipticity > 1: swapping minor <-> major axis")
        ellip = 1/ellip
        theta = np.pi/2 + theta
    sminax = smajax * np.sqrt(1 - ellip)
    apertures = [EllipticalAnnulus(center, a_in=a-a_delta, a_out=a+a_delta, b_out=b, theta=theta) for a, b in zip(smajax, sminax)]
    return apertures

def plot_annuli(data, apertures):
    my_img = plt.imshow(data, origin='lower')
    for ann in apertures:
        ann.plot(color='white')
    plt.colorbar(my_img);
    plt.show()


def compute_stellar_specific_angmom(subsnap, sb_mag, width, resolution, smajax, ellip, a_delta, theta):
    v_los_sim = pynbody.plot.sph.image(subsnap.s, qty='vz', av_z='rho', width=width, resolution=resolution, log=False)
    v_disp_sim = pynbody.plot.sph.image(subsnap.s, qty='v_disp', av_z='rho', width=width, resolution=resolution, log=False)
    lum = to_astropy_quantity(sb_mag, units='mag/arcsec**2')
    v_los = to_astropy_quantity(v_los_sim)
    v_disp = to_astropy_quantity(v_disp_sim)
    lum_annuli = integrate_annulus(lum, center, smajax, ellip, a_delta, theta)
    v_los_annuli = integrate_annulus(v_los, center, smajax, ellip, a_delta, theta)
    v_disp_annuli = integrate_annulus(v_disp, center, smajax, ellip, a_delta, theta)
    stellar_specific_angmom = ss_angmom(lum_annuli, smajax, v_los_annuli, v_disp_annuli)
    return stellar_specific_angmom

logger.info("Opening file")
snap = "/home/michele/sim/MoRIA/M1-10_Verbeke2017/M10sim41001/snapshot_0036"
s = pynbody.load(snap)
max_boxsize = 4000
s.properties['boxsize'] = pynbody.units.Unit("{} kpc".format(max_boxsize))
s.physical_units()
logger.info(s.properties)

width = 5
resolution = 500

logger.info("Set width: {:3d} kpc; resolution {:3d}".format(width, resolution) )
pynbody.analysis.halo.center(s.g)#, vel=False)
subsnap = s[pynbody.filt.Cuboid('{} kpc'.format(-width*1.1))]

logger.info("Computing overall angular momentum")

pynbody.analysis.angmom.faceon(subsnap.s, disk_size=subsnap.s['pos'].max())

logger.info("Computing surface brightness")

fig, ax = plt.subplots()
sb = surface_brightness(subsnap.s, width=width, resolution=resolution, lum_pc2=True, subplot=ax, noplot=not(SHOW))
if SHOW:
    plot_angmom(subsnap.s, ax)
    plt.show()

r_eff_kpc3d = pynbody.analysis.luminosity.half_light_r(subsnap, cylindrical=False)

r_eff_kpc = pynbody.analysis.luminosity.half_light_r(subsnap, cylindrical=True)

r_eff_pix = kpc2pix(r_eff_kpc, width=width, resolution=resolution)

logger.info("Computed R_eff:\n 2D: {:.4f} kpc\n 3D: {:.4f} kpc".format(r_eff_kpc, r_eff_kpc3d))

if FIT_PROFILE:
    logger
    r_bins, sbp = sb_profile(band='v')
    sersic1D = fit_sersic_1D(r_eff_kpc3d, r_bins, sbp, show=SHOW)

    print(sersic1D)

    N_0 = sersic1D.n

logger.info("Fitting Sersic2D")

sersic2D = fit_sersic_2D(sb, r_eff=r_eff_pix, n=N_0, resolution=resolution, ellip=ELLIP_0, theta=THETA_0)

aper = plot_aperture_geometry(sb, sersic=sersic2D, resolution=resolution, show=SHOW)
print(aper)


# Do the photometry

sb_mag = surface_brightness(subsnap.s, width=width, resolution=resolution, lum_pc2=False, noplot=not(SHOW))

center = (sersic2D.x_0.value, sersic2D.y_0.value)
a_delta = 20
smajax = np.arange(30, 200, a_delta)
print(smajax)

logger.info("Created apertures")

apertures = create_apertures(center, smajax, sersic2D.ellip.value, a_delta, sersic2D.theta.value)
if SHOW:
    plot_annuli(sb_mag, apertures)

logger.info("Doing aperture photometry")

sb_ap = to_astropy_quantity(sb_mag)
# flux_table = aperture_photometry(sb_ap, apertures)
# for col in flux_table.colnames:
#     flux_table[col].info.format = '%.8g'  # for consistent table output
# print(flux_table)


p = integrate_annulus(sb_ap, center, smajax, sersic2D.ellip.value, a_delta, sersic2D.theta.value)
print(p)
### 

logger.info("Computing lambda_R")

lambda_R = compute_stellar_specific_angmom(subsnap, sb_mag, width, resolution, smajax, sersic2D.ellip.value, a_delta, sersic2D.theta.value)
print(lambda_R)

###### 
