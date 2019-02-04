import os
import simulation
from dump_features import dump_features
import pickle
import matplotlib.pylab as plt
import numpy as np
import glob

from astropy.table import Table

SIMPATH = '/home/michele/sim/MySimulations/ng'

good_sims = sorted(glob.glob(os.path.join(SIMPATH, 'mb*')))

def get_sim_traj(sim_name):
    s, t = sim_name.split('002_')
    return s+'002', t

if __name__ == '__main__':
    # mega for-loop
    print(good_sims)

    for sim_name in map(os.path.basename, good_sims):
        SIM, TRAJ = get_sim_traj(sim_name)
        sim_path = os.path.join(SIMPATH, "{}_{}".format(SIM, TRAJ), "out")
        print(sim_name)
        sim = simulation.Simulation(sim_path, snap_indexes=slice(None, None, 1))
        load_pickle = False
        radius = 5
        outname = '{}_{}_s{}.fits'.format(SIM, TRAJ, radius)
        if not os.path.isfile(outname):
            dump_features(sim, outname, radius=radius)

    # Read last table
    tbl = Table.read(outname)
        # fig, ax = plt.subplots()
        # im = ax.scatter(mass, sigma_star, c=sim.times)
        # cbar = fig.colorbar(im)
        # ax.set_ylabel('$\sigma$ [km/s]')
        # ax.set_xlabel('$\log_{10} M_\star$  [$M_\odot$]')
        # cbar.ax.set_ylabel('time [Gyr]')
        # m = np.linspace(np.array(mass).min(), np.array(mass).max(), 100)
        # ax.plot(m, no_dm(m), 'r--');
        # ax.set_xscale("log")
        # # ax.set_yscale("log")
        # ax.set_title('{}_{}'.format(SIM, TRAJ));
        # fig.savefig("{}_{}_sig_star.png".format(SIM, TRAJ))
