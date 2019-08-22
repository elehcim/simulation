import os
import simulation
import numpy as np
import pandas as pd
from astropy.table import Table
import matplotlib.pyplot as plt
import pynbody
import tqdm
import argparse
from simulation.util import setup_logger, get_sim_name
import astropy.units  as u
import gc

logger = setup_logger('save_center', logger_level='INFO')

def compute_center(sim_path, vel=False, outname=None):
    sim = simulation.Simulation(sim_path)# , snap_indexes=slice(1, 4))
    times = sim.times
    if outname is None:
        outname = get_sim_name(sim_path) + "_cen.fits"

    nan_arr = np.empty(3, dtype=np.float32)
    nan_arr[:] = np.nan

    cen_list = list()
    vcen_list = list()

    for i, snap in enumerate(tqdm.tqdm(sim)):
        cen_list.append(pynbody.analysis.halo.center(snap.s, retcen=True))
        if vel:
            try:
                pynbody.analysis.halo.center(snap.s, vel=False)
                vcen_new = pynbody.analysis.halo.vel_center(snap.s, retcen=True)
                logger.info("New velocity center: {}".format(vcen_new))
                vcen_list.append(vcen_new)

            except ValueError as e:
                # Usually is 'Insufficient particles around center to get velocity'
                logger.error(e)
                # Get at least the time
                vcen_list.append(nan_arr)
        else:
            vcen_list.append(nan_arr)


    cen = np.array(cen_list)
    vcen = np.array(vcen_list)
    del snap
    sim.snap_list[i] = None  # destroying references to the snap and the list
    gc.collect()

    # write value

    tbl = Table({'t': np.array(times) * u.Gyr,
                 'cx':np.array(cen[:,0]),
                 'cy':np.array(cen[:,1]),
                 'cz':np.array(cen[:,2]),
                 'vcx':np.array(vcen[:,0]),
                 'vcy':np.array(vcen[:,1]),
                 'vcz':np.array(vcen[:,2]),
                 },
                 )
    tbl.write(outname, overwrite=True)
    return sim


def parse_args(cli=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("sim_path", help="Simulation path")
    parser.add_argument("--velocity", '-v', help="Save velocity center", action="store_true")
    parser.add_argument("--outname", "-o", help="Output file name", default=None)
    args = parser.parse_args(cli)
    return args


def main(cli=None):
    args = parse_args(cli)
    compute_center(args.sim_path, args.velocity, args.outname)


if __name__ == '__main__':
    main()


