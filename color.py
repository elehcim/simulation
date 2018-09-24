import pynbody
import numpy as np
import matplotlib.pyplot as plt
from pynbody.plot.stars import convert_to_mag_arcsec2
from pynbody.derived import lum_den_template
from scipy.ndimage.filters import gaussian_filter

def kpc2pix(qty_kpc, width, resolution):
    kpc_per_pixel = width/resolution
    return int(np.floor(qty_kpc/kpc_per_pixel))

def my_convert_to_mag_arcsec2(image):
    assert image.units=="pc^-2"
    pc2_to_sqarcsec = 2.3504430539466191e-09
    img_mag_arcsec2 = -2.5 * np.log10(image * pc2_to_sqarcsec)
    img_mag_arcsec2.units = pynbody.units.arcsec**-2
    return img_mag_arcsec2

def color_plot(snap, bands=('b','i'), width=10, resolution=500, mag_filter=29, subplot=None,
               center=False, title=None, gaussian_sigma=None, cmap_name='seismic', **kwargs):
    """
    Plot the color as defined by the tuple `bands`

    Parameters
    ----------

    gaussian_sigma: in kpc is the sigma of the gaussian to convolve with the image, to make it more realistic

    mag_filter: all region with magnitude/arcsec^2 higher will be set to NaN
    """
    assert len(bands) == 2

    if subplot:
        fig, ax = subplot.figure, subplot
    else:
        fig, ax = plt.gcf(), plt.gca()

    if center:
        pynbody.analysis.halo.center(snap.s, vel=False);

    # create color
    color_name = '{}-{}'.format(*bands)
    snap.s['{}_mag'.format(color_name)] = snap.s['{}_mag'.format(bands[0])] - snap.s['{}_mag'.format(bands[1])]    
    snap.s['{}_lum_den'.format(color_name)] = lum_den_template(color_name, snap.s)

    # plot color in 10^(-0.4) mag per unit surface
    color_pc2 = pynbody.plot.sph.image(snap.s, qty=color_name + '_lum_den', units='pc^-2',
                                       noplot=True, width=width, log=False, resolution=resolution, **kwargs)

    # convert to mag/arcsec**2
    color_mag_arcsec2 = my_convert_to_mag_arcsec2(color_pc2)

    if gaussian_sigma is not None:
        sigma_pix = kpc2pix(gaussian_sigma, width, resolution)
        color_mag_arcsec2 = gaussian_filter(color_mag_arcsec2, sigma_pix)

    # Filter below a certain magnitude
    if mag_filter is not None:
        color_mag_arcsec2[color_mag_arcsec2 > mag_filter] = np.nan
    cmap = plt.get_cmap(cmap_name)
    cmap.set_bad('black')
    img = ax.imshow(color_mag_arcsec2, cmap=cmap, extent=(-width/2, width/2, -width/2, width/2), origin='lower')
    cbar = ax.figure.colorbar(img);
    ax.set_xlabel('x/kpc')
    ax.set_ylabel('y/kpc')
    cbar.set_label('{} [mag/arcsec$^2$]'.format(color_name.upper()));
    # cont = ax.contour(img, cmap='flag', extent=(-width/2, width/2, -width/2, width/2))
    if title is not None:
        ax.set_title(title)
    plt.draw()
    return color_mag_arcsec2
