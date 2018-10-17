import pynbody
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter

def kpc2pix(qty_kpc, width, resolution):
    kpc_per_pixel = width/resolution
    return int(np.floor(qty_kpc/kpc_per_pixel))

def pix2kpc(qty_pix, width, resolution):
    kpc_per_pixel = width/resolution
    return qty_pix*kpc_per_pixel

def convert_to_mag_arcsec2(image):
    """Convert a pynbody.SimArray of luminosity density to one in mag/arcsec^2
    
    At 10 pc (distance for absolute magnitudes), 1 arcsec is 10 AU=1/2.06e4 pc
    In [5]: (np.tan(np.pi/180/3600)*10.0)**2
    Out[5]: 2.3504430539466191e-09
    1 square arcsecond is thus 2.35e-9 pc^2
    """
    pc2_to_sqarcsec = 2.3504430539466191e-09
    img_mag_arcsec2 = -2.5 * np.log10(image.in_units("pc^-2") * pc2_to_sqarcsec)
    img_mag_arcsec2.units = pynbody.units.arcsec**-2
    return img_mag_arcsec2

def surface_brightness(snap, band='v', width=10, resolution=500, mag_filter=29, gaussian_sigma=None,
                       subplot=None, show_cbar=True, center=False, title=None,  cmap_name='viridis', **kwargs):
    """
    Plot and returns the surface brightness in mag/arcsec^2 as defined by `band`.

    It computes the SPH averaged projected luminosity density of the star particles, and
    convert them to mag/arcsec^2

    Parameters
    ----------
    snap : pynbody.SimSnap
        The simulation snapshot
    band : str
        one the available bands: ['u', 'b', 'v', 'r', 'i', 'j', 'h', 'k']
    width : float or pynbody.Unit or str
        if float, the size in kpc of the map. Otherwise use the distance in the given unit
    resolution : int
        number of pixels per side of the image
    center : bool
        if True centers the snap on the star family
    gaussian_sigma : float or None
        in kpc is the sigma of the gaussian to convolve with the image, to make it more realistic
    mag_filter : float or None
        all region with magnitude/arcsec^2 higher will be set to NaN
    kwargs : dict
        optional keyword arguments to be passed to the `pynbody.plot.sph.image` function

    Returns
    -------
    sb_mag_arcsec2 : pynbody.SimArray
        The map of surface brightness in mag/arcsec^2
    """
    if subplot:
        fig, ax = subplot.figure, subplot
    else:
        fig, ax = plt.gcf(), plt.gca()

    if center:
        pynbody.analysis.halo.center(snap.s, vel=False);


    # *_lum_den property is in 10**(-0.4 mag) per unit volume.
    # do a SPH map in 10^(-0.4) mag per unit surface
    pc2 = pynbody.plot.sph.image(snap.s, qty=band + '_lum_den', units='pc^-2',
                                 noplot=True, width=width, log=False, resolution=resolution, **kwargs)

    # convert to mag/arcsec**2
    sb_mag_arcsec2 = convert_to_mag_arcsec2(pc2)

    # Apply the gaussian smoothing
    if gaussian_sigma is not None:
        sigma_pix = kpc2pix(gaussian_sigma, width, resolution)
        sb_mag_arcsec2 = gaussian_filter(sb_mag_arcsec2, sigma_pix)

    # Filter above a certain magnitude
    if mag_filter is not None:
        sb_mag_arcsec2[sb_mag_arcsec2 > mag_filter] = np.nan

    cmap = plt.get_cmap(cmap_name)
    cmap.set_bad('black')

    # Do the plot
    img = ax.imshow(sb_mag_arcsec2, cmap=cmap, extent=(-width/2, width/2, -width/2, width/2), origin='lower')

    if show_cbar:
        cbar = ax.figure.colorbar(img);
        cbar.set_label('{} [mag/arcsec$^2$]'.format(band.upper()));
    ax.set_xlabel('x/kpc')
    ax.set_ylabel('y/kpc')
    # TODO include contour (isophotes), maybe with a paramenter
    # cont = ax.contour(img, cmap='flag', extent=(-width/2, width/2, -width/2, width/2))
    if title is not None:
        ax.set_title(title)
    plt.draw()
    return sb_mag_arcsec2


def color_plot(snap, bands=('b','i'), width=10, resolution=500, mag_filter=29,  gaussian_sigma=None,
               subplot=None, center=False, title=None, cmap_name='seismic',
               vmin=None, vmax=None, **kwargs):
    """
    Plot the color as defined by the tuple `bands`

    It computes the SPH averaged projected luminosity density of the star particles for
    the two bands, convert them to mag/arcsec^2 and subtract them.

    TODO maybe we can use directly the surface_brightness function.

    Parameters
    ----------
    snap : pynbody.SimSnap
        The simulation snapshot
    bands : tuple
        2-tuple of the available bands: ['u', 'b', 'v', 'r', 'i', 'j', 'h', 'k']
    width : float or pynbody.Unit or str
        if float, the size in kpc of the map. Otherwise the distance given
    resolution : int
        number of pixels per side of the image
    center : bool
        if True centers the snap on the star family
    gaussian_sigma : float or None
        in kpc is the sigma of the gaussian to convolve with the image, to make it more realistic
    mag_filter : float or None
        all regions with highr magnitude/arcsec^2 in at least one color map will be set to NaN
    kwargs : dict
        optional keyword arguments to be passed to the `pynbody.plot.sph.image` function

    Returns
    -------
    color_mag_arcsec2 : pynbody.SimArray
        The map of the color in mag/arcsec^2
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
    band0_mag_arcsec2 = convert_to_mag_arcsec2(band0_pc2)
    band1_mag_arcsec2 = convert_to_mag_arcsec2(band1_pc2)

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
        
        # from astropy.convolution import Gaussian2DKernel, convolve
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