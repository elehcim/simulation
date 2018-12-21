import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from snap_io import load_first_last_snap, snapshot_file_list

path = "/home/michele/sim/MySimulations/mb.beta/beta_model.1/out"

def get_times(path):
    snaplist = snapshot_file_list(os.path.expanduser(path), include_dir=True)
    times_map = map(os.path.getmtime, snaplist)
    times = np.fromiter(times_map, dtype=float)
    return times


import gnuplotlib as gp
x = np.arange(len(times))

times = get_times(path) #sys.argv[1])

gp.plot( (x, times), _with  = 'lines', terminal = 'dumb 80,40', unset    = 'grid')