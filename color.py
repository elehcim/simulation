import pynbody
import numpy as np
import matplotlib.pyplot as plt
from astropy.convolution import Gaussian2DKernel, convolve
from pynbody.derived import lum_den_template
from scipy.ndimage.filters import gaussian_filter

def kpc2pix(qty_kpc, width, resolution):
    kpc_per_pixel = width/resolution
    return int(np.floor(qty_kpc/kpc_per_pixel))

def pix2kpc(qty_pix, width, resolution):
    kpc_per_pixel = width/resolution
    return qty_pix*kpc_per_pixel

def my_convert_to_mag_arcsec2(image):
    assert image.units=="pc^-2"
    pc2_to_sqarcsec = 2.3504430539466191e-09  # == 25-5log10(5)
    img_mag_arcsec2 = -2.5 * np.log10(image * pc2_to_sqarcsec)
    img_mag_arcsec2.units = pynbody.units.arcsec**-2
    return img_mag_arcsec2

def color_plot(snap, bands=('b','i'), width=10, resolution=500, mag_filter=29, subplot=None,
               center=False, title=None, gaussian_sigma=None, cmap_name='seismic',
               vmin=None, vmax=None, **kwargs):
    """
    Plot the color as defined by the tuple `bands`

    Parameters
    ----------

    gaussian_sigma: in kpc is the sigma of the gaussian to convolve with the image, to make it more realistic

    mag_filter: all region with magnitude/arcsec^2 higher than mag_filter will be set to NaN
    """
    # create color
    assert len(bands) == 2

    if subplot:
        fig, ax = subplot.figure, subplot
    else:
        fig, ax = plt.gcf(), plt.gca()

    if center:
        pynbody.analysis.halo.center(snap.s, vel=False);

    # plot the two in 10^(-0.4) mag per pc**2
    band0_pc2 = pynbody.plot.sph.image(snap.s, qty=bands[0] + '_lum_den', units='pc^-2', noplot=True, width=width, log=False, resolution=resolution, **kwargs)
    band1_pc2 = pynbody.plot.sph.image(snap.s, qty=bands[1] + '_lum_den', units='pc^-2', noplot=True, width=width, log=False, resolution=resolution, **kwargs)

    # convert to mag/arcsec**2
    band0_mag_arcsec2 = my_convert_to_mag_arcsec2(band0_pc2)
    band1_mag_arcsec2 = my_convert_to_mag_arcsec2(band1_pc2)

    # Filter below a certain magnitude
    if mag_filter is not None:
        band0_mag_arcsec2[band0_mag_arcsec2 > mag_filter] = np.nan
        band1_mag_arcsec2[band1_mag_arcsec2 > mag_filter] = np.nan

    # create color
    color_mag_arcsec2 = band0_mag_arcsec2 - band1_mag_arcsec2
    color_name = '{}-{}'.format(*bands)

    if gaussian_sigma is not None:
        sigma_pix = kpc2pix(gaussian_sigma, width, resolution)
        print("Smoothing with gaussian kernel with sigma = {} pixel".format(sigma_pix))
        color_mag_arcsec2 = gaussian_filter(color_mag_arcsec2, sigma_pix)

        # gaussian_2D_kernel = Gaussian2DKernel(sigma_pix)
        # color_mag_arcsec2 = convolve(color_mag_arcsec2, gaussian_2D_kernel)

    cmap = plt.get_cmap(cmap_name)
    cmap.set_bad('black')
    extent = (-width/2, width/2, -width/2, width/2)
    img = ax.imshow(color_mag_arcsec2, cmap=cmap, interpolation='none', extent=extent, origin='lower', vmin=vmin, vmax=vmax)
    cbar = ax.figure.colorbar(img);
    ax.set_xlabel('x/kpc')
    ax.set_ylabel('y/kpc')
    cbar.set_label('{} [mag/arcsec$^2$]'.format(color_name.upper()));
    # cont = ax.contour(img, cmap='flag', extent=(-width/2, width/2, -width/2, width/2))
    if title is not None:
        ax.set_title(title)
    plt.draw()
    return color_mag_arcsec2