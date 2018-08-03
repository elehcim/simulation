import simulation
import pynbody
import matplotlib.pylab as plt
import numpy as np
import argparse
from pprint import pprint

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

def plot_sfh(sim, ax_sfh, label='SFR', plot_traj=False, force=False, color=None, **kwargs):
    pprint(sim.properties)
    # fig, ax_sfh = plt.subplots(1, figsize=figsize, dpi=dpi)
    # fig.canvas.set_window_title(sim.sim_id)

    sim.plot_sfh(ax_sfh, label=label, **kwargs)

    if plot_traj:
        sim.compute_cog(save_cache=True, force=force)
        ax_r1 = ax_sfh.twinx()
        if sim.is_moving_box:
            ax_r1.plot(sim.trace.t, sim.trace.r, '--', color=color, alpha=0.5, label='$d$')
        else:
            ax_r1.plot(sim.times, np.linalg.norm(sim.cog, axis=0), '--', color=color, alpha=0.5, label='$d$')
        ax_r1.set_ylabel("r [kpc]")
        lines2, labels2 = ax_r1.get_legend_handles_labels()
    else:
        lines2, labels2 = [], []
    # Manage the labels
    lines, labels = ax_sfh.get_legend_handles_labels()
    ax_sfh.legend(lines + lines2, labels + labels2, loc=0)
    # plt.show()

def main(cli=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('sim_dir', nargs='+')
    parser.add_argument('--labels', nargs='+', default=None)
    parser.add_argument('-t', '--traj', action='store_true')
    args = parser.parse_args(cli)

    print(args.sim_dir)
    figsize=(6,4)
    dpi=200
    fig, ax_sfh = plt.subplots(1, figsize=figsize, dpi=dpi)
    if args.labels is None:
        args.labels = ['SFR'] * len(args.sim_dir)
    sims = []
    for sim_dir in args.sim_dir:
        force_cosmo = True if 'MoRIA' in sim_dir else False
        sims.append(simulation.Simulation(sim_dir, force_cosmo=force_cosmo))
    # compute t_max
    t_max = 0
    for sim in sims:
        if sim.times.max() > t_max:
            t_max = sim.times.max()
    for sim, label, color in zip(sims, args.labels, colors):
        print(sim, label, color)
        fig.canvas.set_window_title(sim.sim_id)
        plot_sfh(sim, ax_sfh, label, args.traj, range=[0, t_max], color=color) #, bins=30)
    plt.show()

if __name__ == '__main__':
    main()
