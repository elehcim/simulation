import simulation
import pynbody
import matplotlib.pylab as plt
import numpy as np
import argparse
from pprint import pprint

def plot_sfh(sim, plot_traj=False, figsize=(6,4), dpi=200):
    pprint(sim.properties)
    fig, ax_sfh = plt.subplots(1, figsize=figsize, dpi=dpi)
    fig.canvas.set_window_title(sim.sim_id)

    sim.plot_sfh(ax_sfh, label='SFR')

    if plot_traj:
        sim.compute_cog(save_cache=True, force=True)
        ax_r1 = ax_sfh.twinx()
        ax_r1.plot(sim.times, np.linalg.norm(sim.cog, axis=0), 'r--', alpha=0.5, label='$d$')
        ax_r1.set_ylabel("r [kpc]")
        lines2, labels2 = ax_r1.get_legend_handles_labels()
    else:
        lines2, labels2 = [], []
    # Manage the labels
    lines, labels = ax_sfh.get_legend_handles_labels()
    ax_sfh.legend(lines + lines2, labels + labels2, loc=0)
    plt.show()

def main(cli=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('sim_dir')
    parser.add_argument('-t', '--traj', action='store_true')
    args = parser.parse_args(cli)

    sim = simulation.Simulation(args.sim_dir)
    plot_sfh(sim, args.traj)


if __name__ == '__main__':
    main()