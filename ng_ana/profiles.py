import os
import simulation
import pynbody
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from data_pd.dump_features import dump_features
from simulation.sfh_in_box import binned_sfh, plot_binned_sfh
import tqdm
import astropy.units as u
from simulation.derived import feh, mgfe
from simulation.angmom import sideon

SIM_NUM = 62002

SIM_PATH = "/home/michele/sim/MySimulations/ng/mb.{}_pXX_a800_r600".format(SIM_NUM)
PERI_LIST = [50, 100, 150, 200, 300]

NTH = 1

for peri in PERI_LIST:
    sim_path = os.path.join(SIM_PATH.replace('XX', str(peri)), 'out')
    sim = simulation.Simulation(sim_path, snap_indexes=slice(None, None, NTH))
    surf_bright = {}
    v_circ = {}
    dens = {}
    for i, snap in enumerate(tqdm.tqdm(sim)):
        try:
            # faceon(snap.s)
            sideon(snap.s)
            # print(snap.s)
            snap.properties['eps'] = 0.03
            p = pynbody.analysis.profile.Profile(snap.s, max=4, nbins=100, type='lin', ndim=2)
            surf_bright[i+1] = p['sb']
            v_circ[i+1] = p['v_circ']
            dens[i+1] = p['density']
        except ValueError as e:
            print(e)

    tbl_v_circ = Table({str(k):v*u.km/u.s for (k,v) in v_circ.items()})
    tbl_sb = Table({str(k):v*u.mag/u.arcsec**2 for (k,v) in surf_bright.items()})
    tbl_dens = Table({str(k):v.in_units('Msol kpc**-2')*u.solMass/u.kpc**2 for (k,v) in dens.items()})

    tbl_v_circ['rbins'] = tbl_sb['rbins'] = tbl_dens['rbins'] = p['rbins']*u.kpc

    tbl_v_circ.write('{}p{}prof_v_circ.fits'.format(str(SIM_NUM)[:2], sim.peri), overwrite=True)
    tbl_sb.write('{}p{}prof_sb.fits'.format(str(SIM_NUM)[:2], sim.peri), overwrite=True)
    tbl_dens.write('{}p{}prof_dens.fits'.format(str(SIM_NUM)[:2], sim.peri), overwrite=True)