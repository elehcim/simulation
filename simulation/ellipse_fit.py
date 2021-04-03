import numpy as np
from skimage.measure import find_contours, EllipseModel
import matplotlib.pyplot as plt
import pynbody
from simulation import derived
from collections import namedtuple

ELL_KEYS = ('xc', 'yc', 'a', 'b', 'theta')
EllParams = namedtuple('EllParams', ELL_KEYS)

def ellipsemodel2patch(e, edgecolor='red', **kwargs):
    if edgecolor in kwargs:
        del kwargs['edgecolor']
    return ellipseparams2patch(e.params, edgecolor=edgecolor)


def ellipseparams2patch(params, edgecolor='red', **kwargs):
    if edgecolor in kwargs:
        del kwargs['edgecolor']

    import matplotlib.patches as mpatches
    center = params[0:2]
    a, b = params[2:4]
    theta = params[4]
    return mpatches.Ellipse(center, 2*a, 2*b, theta*180/np.pi, edgecolor=edgecolor, facecolor='none', **kwargs)


def transform_contour(contour, width, resolution):
    """Contour is just a nx2 np.array of (rows, columns), I should transform to physical coordinates"""
    myc = contour*width/resolution-width/2
    # Note that there is a transposition. Result of contour are in (row, column) format
    y = myc[:,0]
    x = myc[:,1]
    return x, y

def plot_contours(cs, width, resolution):
    for c in cs:
        x, y = transform_contour(c, width, resolution)
        plt.plot(x, y)
    plt.gca().set_aspect('equal')
    plt.xlim(-width/2,width/2)
    plt.ylim(-width/2,width/2)


def get_longest_contour(contours):
    longest_contour = np.array([])
    for segment in contours:
        if len(segment) > len(longest_contour):
            longest_contour = segment
    return longest_contour


def fit_ellipse_to_contour(contour):
    ellipse = EllipseModel()
    ellipse.estimate(contour)
    return ellipse

class NoIsophoteFitError(Exception):
    pass

def fit_contour(img, threshold, width, resolution):
    cs = find_contours(img, threshold)
    longest_contour = get_longest_contour(cs)
    if not longest_contour.size:
        return EllParams(*([np.nan]*len(ELL_KEYS)))

    x, y = transform_contour(longest_contour, width, resolution)
    xy = np.vstack([x,y]).T

    ell = fit_ellipse_to_contour(xy)
    # print(ell.params)
    if ell.params is None:
        raise NoIsophoteFitError(f"Can't fit isophote {threshold}")
    return EllParams(*ell.params)

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
    cs = find_contours(rho, threshold)
    longest_contour = get_longest_contour(cs)

    x, y = transform_contour(longest_contour, width, resolution)
    xy = np.vstack([x,y]).T

    ell = fit_ellipse_to_contour(xy)
    print(ell.params)
    ep = ellipsemodel2patch(ell)
    plt.plot(xy[:,0], xy[:,1])
    # plot_contours(cs, width, resolution)
    plt.gca().add_patch(ep)
    plt.show()