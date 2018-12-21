import os
import matplotlib.pyplot as plt
# import simulation
import numpy as np
import pandas as pd
from ..parsers.parse_trace import parse_trace
import matplotlib.ticker as tck


# TODO use argparse
a = np.loadtxt('ssam.69p2.xy.sph/tl.dat')
# a = np.loadtxt('ssam.69p2/tl.dat')

sim_path = '/home/michele/sim/MySimulations/hi_osc/mb.69002_p200_a800_r600/out'
# sim = simulation.Simulation(sim_path)
trace = parse_trace(os.path.join(sim_path, 'trace.txt'))

fig, ax = plt.subplots()
ax.plot(a[:,0], a[:,1])
# FIXME is it in Gyr really?
ax.set_xlabel('t [Gyr]')
ax.set_ylabel(r'$\lambda_R$')

ax_bis = ax.twinx()
# ax_bis.plot(trace.t, trace.r, '--r', alpha=0.5)
ax_bis.set_ylabel('r [kpc]')


ax_bis.plot(trace.t, trace.v, '--g', alpha=0.5)
ax_bis.set_ylabel('v [km/s]')

fig, ax3 = plt.subplots()
im = ax3.scatter(a[:, 2], a[:,1], c=a[:,0])
ax3.set_xlabel('ellipticity')
ax3.set_ylabel(r'$\lambda_R$')

cbar = fig.colorbar(im)
cbar.set_label('t [Gyr]')

fig, ax4 = plt.subplots()
ax4.plot(a[:,0], a[:,3]%np.pi)  # Restrict form 0 and pi
ax4.set_xlabel('t [Gyr]')
ax4.set_ylabel(r'$\theta$')

# ax4.yaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
# ax4.yaxis.set_major_locator(tck.MultipleLocator(base=0.5))


ax4_bis = ax4.twinx()
ax4_bis.plot(trace.t, trace.r, '--r', alpha=0.5)
ax4_bis.set_ylabel('r [kpc]')

plt.show()
