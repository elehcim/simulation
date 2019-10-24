import numpy as np
import pynbody
from ellipse_fit import find_contours, get_longest_contour, transform_contour, fit_ellipse_to_contour

def get_radius_mass_hi_ellipse(rho, width, resolution, threshold=pynbody.units.Unit("1 Msol pc**-2")):
    """
    Integrate the rho_HI map to get the map.
    Fit an ellipse to the threshold contour level to get the radius (major axis of the ellipse)

    rho map should be in Msol pc**-2.

    Return M_HI and r_HI in Msol and kpc respectively.


    Parameters:
    ----------
    width:
        The overall width and height of the plot.
        If ``width`` is a float or an int, then it is assumed to be in kpc.
        It can also be passed in as a string indicating the units,
        i.e. '10 kpc', in which case it converted to kpc.

    Returns:
    --------
    m_HI and r_HI: (tuple)
        HI mass in Msol and HI size in kpc.

    """
    # From e.g. pynbody.plot.sph.image
    from pynbody import units as _units
    if isinstance(width, str) or issubclass(width.__class__, _units.UnitBase):
        if isinstance(width, str):
            width = _units.Unit(width)
        width = width.in_units('kpc')

    width = float(width)
    extent = (-width/2, width/2, -width/2, width/2)
    resolution_element = width * _units.kpc/resolution

    integrand = rho[rho >= threshold].sum()
    m_hi = (integrand * resolution_element**2).in_units("Msol")

    cs = find_contours(rho, threshold)
    longest_contour = get_longest_contour(cs)
    if not longest_contour.size:
        return m_hi, np.nan

    x, y = transform_contour(longest_contour, width, resolution)
    xy = np.vstack([x,y]).T
    ell = fit_ellipse_to_contour(xy)

    # Get semi-major axis
    r_hi = np.max(ell.params[2:4])

    return m_hi, r_hi

if __name__ == '__main__':
    snap = pynbody.load('/home/michele/sim/MySimulations/ok_new_adhoc_or_not_affected/mb.69002_p200_a800_r600_new_adhoc/out/snapshot_0011')
    pynbody.analysis.halo.center(snap.s, vel=False)

    width = 20
    resolution = 1000

    rho = pynbody.plot.image(snap.g,
                             resolution=resolution,
                             qty='rho_HI',
                             cmap='gray',
                             units="Msol pc**-2",
                             width=width,
                             noplot=True,
                            );

    threshold = pynbody.units.Unit("1 Msol pc**-2")
    m, r = get_radius_mass_hi_ellipse(rho, width, resolution, threshold=threshold)
    print(m,r)