import os
import simulation
import pynbody
import numpy as np
import pandas as pd
import tqdm
import pickle
from simulation.sfh_in_box import sfh, plot_sfh, plot_binned_sfh
import glob


gCosmoOrb_rdwarf = 50

def sigma(vel):
    return np.sqrt(((vel - vel.mean(axis=0))**2).mean())


def dump_features(sim, output_name, radius=5, do_sfr=True):
    """
    Write the result of the standard analysis on a pickle file.

    times [Gyr], mass_star [Msol], sigma_star [km/s], sigma_gas [km/s], r_eff [kpc], sfr [Msol/yr]

    """

    sphere = pynbody.filt.Sphere(radius * pynbody.units.kpc)
    mass = list()
    r_eff = list()
    sigma_star = list()
    sigma_gas = list()

    for snap in tqdm.tqdm(sim.snap_list):
        pynbody.analysis.halo.center(snap.s, vel=False)
        mass.append(snap.s[sphere]['mass'].sum().in_units('Msol'))
        sphere_r_eff = pynbody.filt.Sphere(gCosmoOrb_rdwarf * pynbody.units.kpc)
        r_eff.append(pynbody.analysis.luminosity.half_light_r(snap.s[sphere_r_eff]))
        sigma_star.append(sigma(snap.s[sphere]['vel']))
        sigma_gas.append(sigma(snap.g[sphere]['vel']))

    # DUMP
    times = sim.times
    if do_sfr:
        _, sfr = sfh(sim)
    else:
        sfr = np.ones(len(times) + 1) * np.nan

    data = [np.array(times), np.array(mass), np.array(sigma_star),
            np.array(sigma_gas), np.array(r_eff), np.array(sfr)]
    with open(output_name, 'wb') as f:
        pickle.dump(data, f)

def make_df(pickle_file):
    # Gyr, Msol,   km/s    ,  km/s    , kpc  ,Msol/yr
    times, mass_star, sigma_star, sigma_gas, r_eff,   sfr    = pickle.load(open(pickle_file, 'rb'))
    df = pd.DataFrame({'t': times,
                       'mass_star':mass_star,
                       'sigma_star':sigma_star,
                       'sigma_gas':sigma_gas,
                       'r_eff':r_eff,
                       'sfr':sfr})
    return df

if __name__=='__main__':
    moria_list = [62002, 71002, 69002]
    radius = 5
    for moria_name in moria_list:
        moria = simulation.MoriaSim(str(moria_name))
        dump_features(moria, "moria{}.pickle".format(moria_name), radius=radius, do_sfr=False)
    # sim_list = ['mb.62002_p300_a800_r600']
    # SIM_PATH ='/media/michele/My Book/Michele/MySimulations/MovingBox/np'
    # for sim_name in sim_list:
    #     sim = simulation.Simulation(os.path.join(SIM_PATH, sim_name, 'out'))
    #     dump_features(sim, "{}.pickle".format(sim_name))