import os
import simulation
from simulation.util import get_sim_name
import pynbody
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
import tqdm
import argparse
import astropy.units as u

from simulation.parsers.parse_info import parse_info

def Omega(v, a):
    return np.cross(v, a)/((np.linalg.norm(v, axis=1)**2)[:, np.newaxis])


def get_omega_box(sim):
    """ Compute angular velocity of the box

    Returns:
    -------
    om : np.ndarray
        The box angular velocity in units of 1/gadget_time_units = (km/s) / kpc

    """

    dt = parse_info(os.path.join(sim._sim_dir,"info.txt")).dt.values

    v_all = np.vstack([sim.trace.vx, sim.trace.vy, sim.trace.vz]).T

    # Finite difference
    num = v_all[1:]-v_all[0:-1]
    # den = (t[1:]-t[0:-1])[:, np.newaxis]
    # the first difference correspond to idx=2 (third) dt
    den = dt[2:, np.newaxis]
    a_all = num/den

    locations_a = np.digitize(sim.times_header, sim.trace.t[1:], right=True)
    v = v_all[locations_a]
    a = a_all[locations_a]
    om = Omega(v, a)
    return om


def write_omega_box(omega, out_name):
    tbl = Table({'omega': omega * u.km/u.s/u.kpc})
    tbl.write(out_name, overwrite=True)


def parse_args(cli=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(dest='sim_path', help="Path to the simulation snapshots")
    args = parser.parse_args(cli)
    return args


def main(cli=None):
    args = parse_args(cli)
    sim_name = get_sim_name(args.sim_path)
    sim = simulation.Simulation(args.sim_path)
    omega = get_omega_box(sim)
    out_name = sim_name + "_omega.fits"
    print('Writing {}...'.format(out_name))
    write_omega_box(omega, out_name)


if __name__=='__main__':
    main()
