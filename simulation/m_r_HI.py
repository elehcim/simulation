import numpy as np
import tqdm
import matplotlib.pyplot as plt
import pynbody
from collections import defaultdict
from .ellipse_fit import find_contours, get_longest_contour, transform_contour, fit_ellipse_to_contour
from.simdata import get_maps_HI


ELL_KEYS = ('xc', 'yc', 'a', 'b', 'theta')


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

    # FIXME use the astropy quantity units if rho is
    if not isinstance(rho, pynbody.array.SimArray):
        rho = pynbody.array.SimArray(rho, units=pynbody.units.Unit("Msol pc**-2"))

    extent = (-width/2, width/2, -width/2, width/2)
    resolution_element = width * _units.kpc/resolution

    integrand = rho[rho >= threshold].sum()
    m_hi = (integrand * resolution_element**2).in_units("Msol")

    cs = find_contours(rho, threshold)
    longest_contour = get_longest_contour(cs)
    if not longest_contour.size:
        return dict(m_hi=m_hi,
                    r_hi=np.nan,
                    ell=dict(zip(ELL_KEYS, [np.nan]*len(ELL_KEYS))),
                    )

    x, y = transform_contour(longest_contour, width, resolution)
    xy = np.vstack([x,y]).T
    ell = fit_ellipse_to_contour(xy)

    # Get semi-major axis
    r_hi = np.max(ell.params[2:4])

    ell_corr = ell.params.copy()
    # swap axis if needed
    if ell.params[2] < ell.params[3]:
        ell_corr[3], ell_corr[2] = ell.params[2:4]
        ell_corr[4] = ell.params[4] - np.pi/2
    ell_corr[4] = (ell_corr[4] + np.pi/2) % (np.pi) - np.pi/2

    return dict(m_hi=m_hi,
                r_hi=r_hi,
                ell=dict(zip(ELL_KEYS, ell_corr)),
                )


def get_HI_size_mass(sim, width=20, resolution=1000):
    m_list = list()
    r_list = list()
    for snap in tqdm.tqdm(sim):
        pynbody.analysis.halo.center(snap.s, vel=False)
        rho = pynbody.plot.image(snap.g,
                             resolution=resolution,
                             qty='rho_HI',
                             cmap='gray',
                             units="Msol pc**-2",
                             width=width,
                             noplot=True,
                            );
        hi_data = get_radius_mass_hi_ellipse(rho, width, resolution)
        m_list.append(hi_data['m_hi'])
        r_list.append(hi_data['r_hi'])
    m_arr = np.array(m_list)
    r_arr = np.array(r_list)
    return m_arr, r_arr


def get_HI_size_mass_from_maps(maps):
    res = defaultdict(list)
    width = maps.meta['WIDTH']
    resolution = maps.meta['RESOL']
    for rho in tqdm.tqdm(maps['sigma_hi']):
        hi_data = get_radius_mass_hi_ellipse(rho, width, resolution)
        for k in hi_data.keys():
            if k=='ell':
                continue
            res[k].append(hi_data[k])
        for k, v in hi_data['ell'].items():
            res[k].append(v)
    return res

def plot_HI_size_mass(m_arr, r_arr, times, label=None):
    fig, ax = plt.subplots()
    sc = ax.scatter(np.log10(m_arr), np.log10(r_arr),
                c=times,
                s=10,
                label='sim' if label is None else label,
    #             c=range(len(sim_interval.times))
               )
    fig.colorbar(sc);
    m_x = np.linspace(m_arr.min(), m_arr.max(), 100)
    # m_x = np.linspace(1e8, 1e9, 100)
    # ax.scatter(np.log10(m_x), wang16_log10(m_x, mu=0.5, nu=3.54), c='k', s=2);
    ax.plot(np.log10(m_x), wang16_log10(m_x, mu=0.5, nu=3.54));
    rm, rp = wang16_scatter_log10(m_x, mu=0.5, nu=3.54)
    ax.fill_between(np.log10(m_x), rm, rp, alpha=0.4, label='Wang16_mu0.5')
    # ax.set_xscale('log')
    # ax.set_yscale('log');
    ax.set_xlabel('$\log_{10}(m_\mathrm{HI}/$M$_\odot$)');
    ax.set_ylabel('$\log_{10}(r_\mathrm{HI}/$kpc)');
    ax.legend();


def wang16(m):
    mu = 0.506
    nu = 3.293
#     nu = 3.4
    return 10**(mu * np.log10(m) - nu)

def wang16_scatter(m):
    mu = 0.506
    nu = 3.293 + np.log10(2) # D -> r
#     nu = 3.4
    sigma = 0.06
    rm = 10**(mu * np.log10(m) - nu - sigma)
    rp = 10**(mu * np.log10(m) - nu + sigma)
    return rm, rp

def wang16_log10(m, mu=0.506, nu=3.293+np.log10(2)):
    return mu * np.log10(m) - nu

def wang16_scatter_log10(m, mu=0.506, nu=3.293+np.log10(2)):
    sigma = 0.06
    rm = mu * np.log10(m) - nu - sigma
    rp = mu * np.log10(m) - nu + sigma
    return rm, rp


def Begun08_log10(r):
    return 1.96 * np.log10(2 * r) + 6.37


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