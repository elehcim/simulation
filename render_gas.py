import os
import pynbody
from pynbody.plot import sph
import matplotlib.pylab as plt

from notebooks.snap_io import load_moria_sim_and_kicked, load_moria, load_kicked, load_sim

SIMNUMBER = "69002_p200.0_a600.0_r600.0_c8.15"
kicked = True

faceon = True

resolution = 500
width = 20 # kpc
vmin=2e-1
vmax=5e-4;
figsize = (15,15)


def gas_image(sim, **kwargs):
    sim.properties.pop('boxsize', None)
    sim.g['smooth'] /= 2
    pynbody.analysis.halo.center(sim, shrink_factor=0.9)
    try:
        img = sph.image(sim.g, qty="rho", units="g cm^-2", **kwargs)
    finally:
        sim.g['smooth'] *= 2
    return img


def main():
    snap_list = load_kicked(SIMNUMBER) if kicked else load_moria(SIMNUMBER)
    
    folder = "pngs_{}".format(SIMNUMBER)
    os.makedirs(folder, exist_ok=True )
    
    for sim in snap_list:
        snap = int(sim.filename[-4:])
        fig, ax = plt.subplots(figsize=figsize)
        filename = os.path.join(folder,"gas_image_{}_{:03d}.png".format(SIMNUMBER,snap))
        title = '$t={:5.2f}$ Gyr, snap={:03d}'.format(sim.properties['time'].in_units("Gyr"), snap)
        gas_image(sim, width=width, vmin=vmin, vmax=vmax, title=title, filename=filename);
        print("Saved", filename)

    
if __name__ == '__main__':
    main()