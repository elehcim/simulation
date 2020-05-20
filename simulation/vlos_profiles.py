# compute maximum LOS velocity from velocity map on the plane of the orbit
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from skimage.measure import profile_line
from simulation.simdata import get_vlos_map


def plot_maps(vlos_map, idx):
    res = vlos_map.shape[1]
    p1, p2 = (0,res//2), (res,res//2)
    plt.imshow(vlos_map[idx], origin='lower')
    line = matplotlib.lines.Line2D(xdata=(p1[0],p2[0]), ydata=(p1[1],p2[1]),linestyle='dashed',color='k')
    plt.axes().add_line(line)
    plt.colorbar()

    # plt.plot(np.abs(prof));


def compute_profile(vlos_map):
    """Compute the profile using the horizontal line on the middle of the image"""

    #TODO do an average on some pixels around center
    res = vlos_map.shape[1]
    prof = profile_line(vlos_map.transpose(1,2,0), (res//2,0), (res//2,res-1))
    return prof


def plot_profile(prof, ax=None):
    """Plot the profile"""
    if ax is None:
        ax = plt
    # ax.plot(np.abs(prof, axis=0))
    ax.plot(prof)


# Entry point function
def get_max_vlos(sim_name, orbit_sideon=True):
    """Return an array of the maximum v_los as computed from the middle horizontal slice of the vlos_map"""
    vlos_map = get_vlos_map(sim_name, orbit_sideon)
    prof = compute_profile(vlos_map)
    max_vlos = np.max(np.abs(prof), axis=0)
    return max_vlos