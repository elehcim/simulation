import pynbody
import numpy as np
import copy
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter

# Values of pynbody (pynbody.luminosity.halo_lum) from http://www.ucolick.org/~cnaw/sun.html which does not exist anymore.
# SUN_ABS_MAGNITUDES = {'u':5.56, 'b':5.45, 'v':4.8, 'r':4.46, 'i':4.1, 'j':3.66, 'h':3.32, 'k':3.28}

# The following values from here: http://mips.as.arizona.edu/~cnaw/sun.html
SUN_ABS_MAGNITUDES = {'u':5.61, 'b':5.44, 'v':4.81,
                      'r':4.43, 'i':4.10,
                      'j':3.67, 'h':3.32, 'k':3.27,
                      'sdss_u':5.49, 'sdss_g':5.23, 'sdss_r':4.53, 'sdss_i':4.19, 'sdss_z':4.01}


def kpc2pix(qty_kpc, width, resolution):
    kpc_per_pixel = width/resolution
    return int(np.floor(qty_kpc/kpc_per_pixel))


def pix2kpc(qty_pix, width, resolution):
    kpc_per_pixel = width/resolution
    return qty_pix*kpc_per_pixel


def convert_to_mag_arcsec2(image, band):
    """Convert a pynbody.SimArray of luminosity density to one in mag/arcsec^2

    At 10 pc (distance for absolute magnitudes), 1 arcsec is 10 AU = 4.848e-5 pc = (1/2.06e4) pc
    In [4]: (np.tan(np.pi/180/3600)*10.0)
    Out[4]: 4.848e-5
    In [5]: (np.tan(np.pi/180/3600)*10.0)**2
    Out[5]: 2.3504430539466191e-09
    1 square arcsecond is thus 2.35e-9 pc^2
    """
    pc2_to_sqarcsec = 2.3504430539466191e-09
    img_mag_arcsec2 = SUN_ABS_MAGNITUDES[band] - 2.5 * np.log10(image.in_units("pc^-2") * pc2_to_sqarcsec)
    img_mag_arcsec2.units = pynbody.units.arcsec**-2
    # equivalently: -2.5*log10(2.3504430539466191e-09) = 21.572
    # so: img_mag_arcsec2 = SUN_ABS_MAGNITUDES[band] + 21.572 - 2.5 * np.log10(image.in_units("pc^-2"))
    # as in wikipedia: https://en.wikipedia.org/wiki/Surface_brightness#Relationship_to_physical_units
    return img_mag_arcsec2


def surface_brightness(snap, band='v', width=10, resolution=500, center=False, lum_pc2=False, noplot=False,
                       mag_filter=None, gaussian_sigma=None, subplot=None, show_cbar=True, cax=None,
                       cmap_name='viridis', title=None, isophotes=0, label_contour=True, bad_pixels='black', **kwargs):
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
    lum_pc2 : bool
        If true return the image in solar luminosity in the given band per pc**2. The default is False
    gaussian_sigma : float or None
        in kpc is the sigma of the gaussian to convolve with the image, to make it more realistic
    mag_filter : float or None
        all region with magnitude/arcsec^2 higher will be set to NaN
    isophotes : int, iterable
        Number of contour levels. It is active only if lum_pc2=False. Default to 0 (no isophotes plot)
    kwargs : dict
        optional keyword arguments to be passed to the `imshow` function

    Returns
    -------
    sb : pynbody.SimArray
        The map of surface brightness in mag/arcsec^2 or in Lsol/pc**2
    """
    if center:
        pynbody.analysis.halo.center(snap.s, vel=False)

    lum_density_name = band + '_lum_density'
    sun_abs_mag = SUN_ABS_MAGNITUDES[band]
    snap.s[lum_density_name] = (10 ** (-0.4 * (snap.s[band + "_mag"] - sun_abs_mag))) * snap.s['rho'] / snap.s['mass']

    snap.s[lum_density_name].units = snap.s['rho'].units/snap.s['mass'].units

    try:
        # snap.s['smooth'] /= 2
        # Projection, units imply projection
        pc2 = pynbody.plot.sph.image(snap.s, qty=lum_density_name, units='pc^-2',
                                 noplot=True, width=width, resolution=resolution)
    finally:
        # snap.s['smooth'] *= 2
        pass

    if lum_pc2:
        sb = pc2
    else:
        sb = convert_to_mag_arcsec2(pc2, band)

    # Apply the gaussian smoothing
    if gaussian_sigma is not None:
        sigma_pix = kpc2pix(gaussian_sigma, width, resolution)
        sb = gaussian_filter(sb, sigma_pix)

    # Filter above a certain magnitude
    if not lum_pc2 and mag_filter is not None:
        sb[sb > mag_filter] = np.nan

    log = kwargs.get('log', False)
    if log and lum_pc2:
        sb = np.log10(sb)

    if noplot:
        return sb

    # Do the plot
    if subplot:
        fig, ax = subplot.figure, subplot
    else:
        fig, ax = plt.gcf(), plt.gca()

    # Copy to remove deprecation warning
    cmap = copy.copy(plt.get_cmap(cmap_name))
    cmap.set_bad(bad_pixels)

    extent = (-width/2, width/2, -width/2, width/2)
    img = ax.imshow(sb, cmap=cmap, extent=extent, origin='lower', **kwargs)


    if show_cbar:
        if lum_pc2:
            cbar_label = '${0}I_{1}$ (L$_{{\odot,{1}}}$/pc$^2$)'.format("Log" if log else "", band.upper())
        else:
            cbar_label = '$\mu_{{{}}}$ (mag/arcsec$^2$)'.format(band.upper())

        from mpl_toolkits.axes_grid1.axes_grid import CbarAxes

        if isinstance(cax, CbarAxes):
            cbar = cax.colorbar(img)
            cbar.set_label_text(cbar_label)
        elif cax:
            cbar = ax.figure.colorbar(img, cax=cax)
            cbar.set_label(cbar_label)
        else:
            cbar = ax.figure.colorbar(img)
            cbar.set_label(cbar_label)

    ax.set_xlabel('x/kpc')
    ax.set_ylabel('y/kpc')
    if not lum_pc2 and isophotes:
        if isinstance(isophotes, int):
            levels = np.linspace(sb.min(), sb.max(), isophotes, dtype=np.int)
        else:
            levels = isophotes
        cont = ax.contour(sb, levels=levels, extent=extent) #  cmap='flag' # very visible countours
        if label_contour:
            ax.clabel(cont, inline=1, fmt='%1.0f')
    if title is not None:
        ax.set_title(title)
    plt.draw()
    return sb


def color_plot(snap, bands=('b','i'), width=10, resolution=500, mag_filter=29, gaussian_sigma=None,
               subplot=None, center=False, title=None, cmap_name='seismic', noplot=False,
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
    color : pynbody.SimArray
        The map of the color in mag/arcsec^2
    """
    # create color
    assert len(bands) == 2

    if center:
        pynbody.analysis.halo.center(snap.s, vel=False)

    sb = list()
    for band in bands:
        sb.append(surface_brightness(snap.s, band=band, width=width, resolution=resolution, lum_pc2=False,
                                     center=False, noplot=True, mag_filter=mag_filter))

    # create color
    color = sb[0] - sb[1]

    if gaussian_sigma is not None:
        sigma_pix = kpc2pix(gaussian_sigma, width, resolution)
        print("Smoothing with gaussian kernel with sigma = {} pixel".format(sigma_pix))
        color = gaussian_filter(color, sigma_pix)

        # from astropy.convolution import Gaussian2DKernel, convolve
        # gaussian_2D_kernel = Gaussian2DKernel(sigma_pix)
        # color = convolve(color, gaussian_2D_kernel)

    if noplot:
        return color

    if subplot:
        fig, ax = subplot.figure, subplot
    else:
        fig, ax = plt.gcf(), plt.gca()

    cmap = plt.get_cmap(cmap_name)
    cmap.set_bad('black')
    extent = (-width/2, width/2, -width/2, width/2)
    img = ax.imshow(color, cmap=cmap, interpolation='none', extent=extent, origin='lower', vmin=vmin, vmax=vmax)
    cbar = ax.figure.colorbar(img)
    ax.set_xlabel('x/kpc')
    ax.set_ylabel('y/kpc')
    color_name = '{}-{}'.format(*bands)
    cbar.set_label('{} (mag/arcsec$^2$)'.format(color_name.upper()))
    # cont = ax.contour(img, cmap='flag', extent=(-width/2, width/2, -width/2, width/2))
    if title is not None:
        ax.set_title(title)
    plt.draw()
    return color
