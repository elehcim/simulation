import sys
import os
import simulation
import matplotlib.pyplot as plt
from parse_trace import parse_trace
from pprint import pprint

MB_BETA_SIM_PATH = '/media/michele/My Book/Michele/MySimulations/MovingBox/beta/'

sim_n = '62002'
peri = 300
traj = 'p{}_a800_r600'.format(peri)

sim_path = os.path.join(MB_BETA_SIM_PATH,'mb.{}'.format(sim_n), traj, 'out')

moria = simulation.MoriaSim(sim_n)
# moria.compute_cog(save_cache=True, force=True)
sim = simulation.Simulation(sim_path)
pprint(moria.properties)
pprint(sim.properties)

figsize = (6,4)
bins=100
fig, ax_sfh = plt.subplots(1, figsize=figsize, dpi=200)
sim.plot_sfh(ax_sfh, label='SFR', bins=bins, last_snap=-1)
moria.plot_sfh(ax_sfh, label="MoRIA", alpha=0.3, bins=bins)

ax_r1 = ax_sfh.twinx()
ax_r1.plot(sim.trace.t, sim.trace.r, 'r--', alpha=0.5, label='$d$')
ax_r1.set_ylabel("r [kpc]")

# Manage the labels
lines, labels = ax_sfh.get_legend_handles_labels()
lines2, labels2 = ax_r1.get_legend_handles_labels()
ax_sfh.legend(lines + lines2, labels + labels2, loc=0)
fig.show()