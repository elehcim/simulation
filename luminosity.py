import pynbody
import numpy as np
import matplotlib.pyplot as plt
from pynbody.plot.stars import convert_to_mag_arcsec2
from pynbody.derived import lum_den_template
from scipy.ndimage.filters import gaussian_filter
from color import kpc2pix, my_convert_to_mag_arcsec2

def luminosity_plot(snap, band='v', width=10, resolution=500, mag_filter=29, subplot=None,
               center=False, title=None, gaussian_sigma=None, cmap_name='viridis', **kwargs):
    """
    Plot the luminosity in mag/arcsec^2 as defined by the band

    Parameters
    ----------

    gaussian_sigma: in kpc is the sigma of the gaussian to convolve with the image, to make it more realistic

    mag_filter: all region with magnitude/arcsec^2 higher will be set to NaN
    """
    if subplot:
        fig, ax = subplot.figure, subplot
    else:
        fig, ax = plt.gcf(), plt.gca()

    if center:
        pynbody.analysis.halo.center(snap.s, vel=False);


    # plot color in 10^(-0.4) mag per unit surface
    pc2 = pynbody.plot.sph.image(snap.s, qty=band + '_lum_den', units='pc^-2',
                                 noplot=True, width=width, log=False, resolution=resolution, **kwargs)

    # convert to mag/arcsec**2
    mag_arcsec2 = my_convert_to_mag_arcsec2(pc2)

    if gaussian_sigma is not None:
        sigma_pix = kpc2pix(gaussian_sigma, width, resolution)
        mag_arcsec2 = gaussian_filter(arcsec2, sigma_pix)

    # Filter below a certain magnitude
    if mag_filter is not None:
        mag_arcsec2[mag_arcsec2 > mag_filter] = np.nan
    cmap = plt.get_cmap(cmap_name)
    cmap.set_bad('black')
    img = ax.imshow(mag_arcsec2, cmap=cmap, extent=(-width/2, width/2, -width/2, width/2), origin='lower')
    cbar = ax.figure.colorbar(img);
    ax.set_xlabel('x/kpc')
    ax.set_ylabel('y/kpc')
    cbar.set_label('{} [mag/arcsec$^2$]'.format(band.upper()));
    # cont = ax.contour(img, cmap='flag', extent=(-width/2, width/2, -width/2, width/2))
    if title is not None:
        ax.set_title(title)
    plt.draw()
    return mag_arcsec2