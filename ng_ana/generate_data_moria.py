import os
import simulation
import pynbody
import pandas as pd
import matplotlib.pylab as plt
import numpy as np
import tqdm
import pickle

gCosmoOrb_rdwarf = 50

def sigma(vel):
    return np.sqrt(((vel - vel.mean(axis=0))**2).mean())

def dump_features(sim, output_name):
    sphere = pynbody.filt.Sphere(5 * pynbody.units.kpc)
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
    data = [times, mass, sigma_star, sigma_gas, r_eff]
    with open(output_name, 'wb') as f:
        pickle.dump(data, f)

if __name__=='__main__':
    # moria_list = [62002, 71002, 69002]
    # for moria_name in moria_list:
    #     moria = simulation.MoriaSim(str(moria_name))
    #     dump_features(moria, "moria{}.pickle".format(moria_name))
    sim_list = ['mb.62002_p300_a800_r600']
    SIM_PATH ='/media/michele/My Book/Michele/MySimulations/MovingBox/np'
    for sim_name in sim_list:
        sim = simulation.Simulation(os.path.join(SIM_PATH, sim_name, 'out'))
        dump_features(sim, "{}.pickle".format(sim_name))