import matplotlib.pylab as plt
import pynbody
from pynbody.plot import sph

def plot_gas(sim, i):
    sim.properties.pop('boxsize', None)
    sim.g['smooth'] /= 2
    pynbody.analysis.halo.center(sim, shrink_factor=0.9)
    fig, ax = plt.subplots(figsize=(10,10))
    sph.image(sim.g, qty="rho", resolution=200, units="g cm^-2",
              width=10, )#, vmin=2e-1, vmax=5e-4,
              # filename="{}/{}_{:03d}".format(foldername, SIMNUMBER, i));
    time = '$t={:5.2f}$ Gyr, snap={:03d}'.format(sim.properties['time'].in_units("Gyr"), i)
    #ax = plt.gca()
    #ax.set_title(time, fontsize=12)
    return fig, ax