import simulation
import pynbody
import matplotlib.pylab as plt
import numpy as np

# sim = simulation.MoriaSim('69002_p200.0_a600.0_r600.0_c8.15_z0', kicked=True)
sim_no_kick = simulation.MoriaSim('69002', kicked=False)

# sim_name = '71002'  # This has the jump in SFR.
sim_name = '69002_p10.0_a600.0_r600.0_c8.15_z0'
sim = simulation.MoriaSim(sim_name, kicked=True)
# sim_no_kick = simulation.MoriaSim(sim_name, kicked=False)




sim.compute_cog(save_cache=True, force=True)


figsize = (6,4)
fig, ax_sfh = plt.subplots(1, figsize=figsize, dpi=200)
sim.plot_sfh(ax_sfh, label='SFR')
sim_no_kick.plot_sfh(ax_sfh, label="MoRIA", alpha=0.3)

ax_r1 = ax_sfh.twinx()
ax_r1.plot(sim.times, np.linalg.norm(sim.cog, axis=0), 'r--', alpha=0.5, label='$d$')
ax_r1.set_ylabel("r [kpc]")

# Manage the labels
lines, labels = ax_sfh.get_legend_handles_labels()
lines2, labels2 = ax_r1.get_legend_handles_labels()
ax_sfh.legend(lines + lines2, labels + labels2, loc=0)
fig.show()