import os
import simulation
import pynbody
import pandas as pd
import matplotlib.pylab as plt
import numpy as np
import ipywidgets
import tqdm
import pickle
from simulation.sfh_in_box import sfh, plot_sfh, plot_binned_sfh
import glob
from astropy import constants as c

SIMPATH = '/home/michele/sim/MySimulations/ng'
DATA_PATH = '.'

G = c.G.to('kpc solMass**-1 km**2 s**-2')

def sigma_no_dm(m_star, r_eff):
    """Compute sigma from the virial radius assuming no dm.
    m_star in 1e10Msol
    r_eff in kpc
    """
    s = 0.45 * (G.value * m_star / r_eff)
    return s


if __name__ == '__main__':
    # mega for-loop
    for sim_name in map(os.path.basename, good_sims):
        SIM, TRAJ = get_sim_traj(sim_name)
        sim_path = os.path.join(SIMPATH, "{}_{}".format(SIM, TRAJ), "out")
        print(sim_name)
        sim = simulation.Simulation(sim_path, snap_indexes=slice(None, None, 1))
        load_pickle = True
        radius = 5
        outname = os.path.join(DATA_PATH, '{}_{}_s{}.pickle'.format(SIM, TRAJ, radius))
        if load_pickle:
            times, mass, sigma_star, sigma_gas, r_eff, sfr = pickle.load(open(outname, 'rb'))
        fig, ax = plt.subplots()
        im = ax.scatter(mass, sigma_star, c=sim.times)
        cbar = fig.colorbar(im)
        ax.set_ylabel('$\sigma$ [km/s]')
        ax.set_xlabel('$\log_{10} M_\star$  [$M_\odot$]')
        cbar.ax.set_ylabel('time [Gyr]')
        m = np.linspace(np.array(mass).min(), np.array(mass).max(), 100)
        ax.plot(m, no_dm(m), 'r--');
        ax.set_xscale("log")
        # ax.set_yscale("log")
        ax.set_title('{}_{}'.format(SIM, TRAJ));
        fig.savefig("{}_{}_sig_star.png".format(SIM, TRAJ))