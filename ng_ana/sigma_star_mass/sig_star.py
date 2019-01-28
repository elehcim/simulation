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
from dump_features import dump_features, sigma
import glob

def no_dm(m):
    s = 5e-4 * (m ** 0.5)
    return s

SIMPATH = '/media/michele/My Book/Michele/MySimulations/MovingBox/np'

good_sims = glob.glob(os.path.join(SIMPATH, 'mb*'))
print(good_sims)
def get_sim_traj(sim_name):
    s, t = sim_name.split('002_')
    return s+'002', t

#SIM = 'mb.71002'
#TRAJ = 'p100_a800_r600'


for sim_name in map(os.path.basename, good_sims):
    SIM, TRAJ = get_sim_traj(sim_name)
    sim_path = os.path.join(SIMPATH, "{}_{}".format(SIM, TRAJ), "out")
    print(sim_name)
    sim = simulation.Simulation(sim_path, snap_indexes=slice(None, None, 1))
    load_pickle = False
    radius = 5
    outname = '{}_{}_s{}.pickle'.format(SIM, TRAJ, radius)
    if load_pickle:
        times, mass, sigma_star, sigma_gas, r_eff = pickle.load(open(outname, 'rb'))
    else:
        if not os.path.isfile(outname):
            dump_features(sim, outname, radius=radius)
        times, mass, sigma_star, sigma_gas, r_eff = pickle.load(open(outname, 'rb'))

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
