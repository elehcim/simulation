import os
import simulation
import pynbody
import numpy as np
import pandas as pd
import tqdm
import pickle
from simulation.sfh_in_box import sfh, plot_sfh, plot_binned_sfh
import glob
from astropy.table import Table
import astropy.units as u

# gCosmoOrb_rdwarf = 50  # kpc
# R_EFF_BORDER = gCosmoOrb_rdwarf
R_EFF_BORDER = 10  # Same as hyplot

def sigma(vel):
    return np.sqrt(((vel - vel.mean(axis=0))**2).mean())


def dump_features(sim, output_name, radius=R_EFF_BORDER, do_sfr=True):
    """
    Write the result of the standard analysis on a pickle file.

    times [Gyr], mass_star [Msol], sigma_star [km/s], sigma_gas [km/s], r_eff [kpc], sfr [Msol/yr]

    """

    sphere = pynbody.filt.Sphere(radius * pynbody.units.kpc)
    mass_star = list()
    r_eff = list()
    r_eff3d = list()
    sigma_star = list()
    sigma_gas = list()
    metals_star = list()

    for snap in tqdm.tqdm(sim.snap_list):
        pynbody.analysis.halo.center(snap.s, vel=False)
        star_m = snap.s[sphere]['mass'].sum().in_units('Msol')
        mass_star.append(star_m)
        r_eff.append(pynbody.analysis.luminosity.half_light_r(snap.s[sphere],  cylindrical=True))
        r_eff3d.append(pynbody.analysis.luminosity.half_light_r(snap.s[sphere], cylindrical=False))
        sigma_star.append(sigma(snap.s[sphere]['vel']))
        sigma_gas.append(sigma(snap.g[sphere]['vel']))
        metals_star.append((snap.s[sphere]['metals'] * snap.s[sphere]['mass']).sum()/star_m)
        #TODO

    # DUMP
    times = sim.times
    if do_sfr:
        _, sfr = sfh(sim)
    else:
        sfr = np.ones(len(times) + 1) * np.nan


    tbl = Table({'t': np.array(times) * u.Gyr,
                 'mass_star':np.array(mass_star) * u.solMass,
                 'sigma_star':np.array(sigma_star) * u.km/u.s,
                 'sigma_gas':np.array(sigma_gas) * u.km/u.s,
                 'r_eff':np.array(r_eff) * u.kpc,
                 'r_eff3d':np.array(r_eff3d) * u.kpc,
                 'sfr':np.array(sfr) * u.solMass/u.yr,
                 'metals_star':np.array(metals_star),
                 },
                 )
    tbl.write(output_name, format='fits', overwrite=True)

# def make_df(pickle_file):
#     # Gyr, Msol,   km/s    ,  km/s    , kpc  ,Msol/yr
#     times, mass_star, sigma_star, sigma_gas, r_eff,   sfr    = pickle.load(open(pickle_file, 'rb'))
#     df = pd.DataFrame({'t': times,
#                        'mass_star':mass_star,
#                        'sigma_star':sigma_star,
#                        'sigma_gas':sigma_gas,
#                        'r_eff':r_eff,
#                        'sfr':sfr})
#     return df

# def make_astropy_table(pickle_file):
#     times, mass_star, sigma_star, sigma_gas, r_eff = pickle.load(open(pickle_file, 'rb'))
#     tbl = Table({'t': times * u.Gyr,
#                  'mass_star':mass_star * u.solMass,
#                  'sigma_star':sigma_star * u.km/u.s,
#                  'sigma_gas':sigma_gas * u.km/u.s,
#                  'r_eff':r_eff * u.kpc,
#                  # 'sfr':sfr * u.solMass/u.yr,
#                  })
#     return tbl

if __name__=='__main__':
    moria_list = [62002, 71002, 69002]
    for moria_name in moria_list:
        moria = simulation.MoriaSim(str(moria_name))
        dump_features(moria, "moria{}.fits".format(moria_name), do_sfr=False)
    # sim_list = ['mb.62002_p300_a800_r600']
    # SIM_PATH ='/media/michele/My Book/Michele/MySimulations/MovingBox/np'
    # for sim_name in sim_list:
    #     sim = simulation.Simulation(os.path.join(SIM_PATH, sim_name, 'out'))
    #     dump_features(sim, "{}.pickle".format(sim_name))