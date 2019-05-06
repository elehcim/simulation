import os
import simulation
import pynbody
import pandas as pd
import matplotlib.pylab as plt
import matplotlib
import numpy as np
import ipywidgets
import tqdm
from simulation.metals import *
from astropy import constants as c
from astropy.table import Table
from mpl_toolkits.axes_grid1 import make_axes_locatable

# SIM_PATH = "/media/michele/My Book/Michele/MySimulations/MovingBox/np/mb.62002_pXX_a800_r600"
SIM_PATH = "/home/michele/sim/MySimulations/ng/mb.71002_pXX_a800_r600"
PERI_LIST = [50, 100, 150, 200, 300]
NTH = 5
sims = dict()
for peri in PERI_LIST:
    sim_path = os.path.join(SIM_PATH.replace('XX', str(peri)), 'out')
    sims["71p{}".format(peri)] = simulation.Simulation(sim_path, snap_indexes=slice(None, None, NTH))

sims['moria'] = simulation.MoriaSim(71002, snap_indexes=slice(65, None, NTH))

fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(12,10))

for ax, sim, title in zip(axs.flatten(), sims.values(), sims.keys()):
    mappable = plt.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=sim.times_header.min(), vmax=sim.times_header.max()))
    for snap in sim:
        plot_mdf(snap.s, ax=ax, nbins=20, color=mappable.to_rgba(snap.header.time), alpha=0.4)
    mappable.set_array([])
    ax.set_ylim(0, 1)
    ax.set_title(title)
    ax.set_xlabel('[Fe/H]')
    ax.set_ylabel('PDF')
    ax.grid()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(mappable, cax=cax, orientation='vertical');