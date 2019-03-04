import os
import simulation
import pynbody
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
import tqdm

from simulation.parsers.parse_info import parse_info

def Omega(v, a):
    return np.cross(v, a)/((np.linalg.norm(v, axis=1)**2)[:, np.newaxis])

def get_omega_box(sim):
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
    Om = Omega(v, a)
    return Om