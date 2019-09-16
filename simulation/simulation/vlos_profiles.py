# compute maximum LOS velocity from velocity map on the plane of the orbit
import os
import matplotlib
import simulation
import matplotlib.pyplot as plt
import pynbody
import numpy as np
import pandas

from astropy.table import Table
from skimage.measure import profile_line
from pynbody import units
from pynbody.filt import Filter

data_path = '/home/michele/sim/analysis/ng_ana/data/'

class Slit(Filter):
    """ Emulate a Slit, so taking a slice at z=const of a snap"""
    def __init__(self, x1, y1=None, x2=None, y2=None):
        self._descriptor = "slit"
        x1, y1, x2, y2 = [units.Unit(x) if isinstance(x, str) else x for x in (x1, y1, x2, y2)]
        if y1 is None:
            y1 = x1
        if x2 is None:
            x2 = -x1
        if y2 is None:
            y2 = -y1

        self.x1, self.y1, self.x2, self.y2 = x1, y1, x2, y2

    def __call__(self, sim):
        x1, y1, x2, y2 = [x.in_units(sim["pos"].units, **sim["pos"].conversion_context())
                                  if units.is_unit_like(x) else x
                                  for x in (self.x1, self.y1, self.x2, self.y2)]

        return ((sim["x"] > x1) * (sim["x"] < x2) * (sim["y"] > y1) * (sim["y"] < y2))

    def __repr__(self):
        x1, y1, x2, y2 = ["'%s'" % str(x)
                                  if units.is_unit_like(x) else x
                                  for x in (self.x1, self.y1, self.x2, self.y2)]
        return "Slit(%s, %s, %s, %s)" % (x1, y1, x2, y2)



# sim_name = simulation.util.get_sim_name(simpath)


def get_vlos_map(sim_name, orbit_sideon):
    if orbit_sideon:
        tbl = Table.read(os.path.join(data_path, 'maps_orbit_sideon_sig_los', sim_name+'_maps.fits'))
    else:
        tbl = Table.read(os.path.join(data_path, 'maps_orbit_faceon_sig_los', sim_name+'_maps.fits'))
    vlos_map = tbl['vlos']
    return vlos_map

def plot_maps(vlos_map, idx):
    res = vlos_map.shape[1]
    p1, p2 = (0,res//2), (res,res//2)
    plt.imshow(vlos_map[idx], origin='lower')
    line = matplotlib.lines.Line2D(xdata=(p1[0],p2[0]), ydata=(p1[1],p2[1]),linestyle='dashed',color='k')
    plt.axes().add_line(line)
    plt.colorbar()

    # plt.plot(np.abs(prof));


def compute_profile(vlos_map):
    #TODO do an average on some pixels around center
    res = vlos_map.shape[1]
    prof = profile_line(vlos_map.transpose(1,2,0), (res//2,0), (res//2,res-1))
    return prof


def plot_profile(prof, ax=None):
    if ax is None:
        ax = plt
    ax.plot(np.max(np.abs(prof), axis=0))

# Entry point function
def get_max_vlos(sim_name, orbit_sideon=True):
    """Return an array of the maximum v_los as computed from the middle horizontal slice of the vlos_map"""
    vlos_map = get_vlos_map(sim_name, orbit_sideon)
    prof = compute_profile(vlos_map)
    max_vlos = np.max(np.abs(prof), axis=0)
    return max_vlos