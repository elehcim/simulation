import pynbody
import numpy as np
import os
from scipy.io import savemat
from simulation import derived
import argparse

def convert_snapshot_to_mat(snap_name, density_threshold, outfile_name):
    """Take a snapshot, cut on the density and provide a mat file with particle positions (x,y,z) and IDs"""
    snap = pynbody.load(snap_name, ignore_cosmo=True)
    gas = snap.gas
    filt = pynbody.filt.HighPass('rho', density_threshold)

    pos = gas[filt]['pos'].view(np.ndarray)
    ids = gas[filt]['iord'].view(np.ndarray)
    rho = gas[filt]['rho'].view(np.ndarray)

    data = dict(x=pos[:,0], y=pos[:,1], z=pos[:,2], id=ids, rho=rho)
    print('output file', outfile_name)
    savemat(outfile_name, data)
    return pos

QUANTITIES = ('temp', 'u',
              'rho', 'mass', 'p',
              'feh', 'mgfe', 'zsph',
              'mass_HI', 'neutral_fraction',
              'vx', 'vy', 'vz', 'v_norm', 'mach', 'cii', 'smooth',
              )

DENSITY_THRESHOLD = 6e-6

def convert_to_mat_all_info(snap_name, density_threshold, outfile_name, quantities=QUANTITIES):
    snap = pynbody.load(snap_name, ignore_cosmo=True)
    gas = snap.gas
    filt = pynbody.filt.HighPass('rho', density_threshold)
    data = dict()
    for qty in list(quantities) + ['iord', 'x', 'y', 'z']:

       data[qty] = gas[filt][qty].view(np.ndarray)
       # dict(x=pos[:,0], y=pos[:,1], z=pos[:,2], id=ids, rho=rho)
    print('output file', outfile_name)
    # print(data)
    savemat(outfile_name, data)

def parse_args(cli=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(dest='snap_name', help="Path to the simulation snapshot")
    parser.add_argument('-t', '--density-threshold', help='Density threshold', type=float, default=DENSITY_THRESHOLD)
    parser.add_argument('-o', '--output', dest='outfile_name', help='Output filename', default=None)
    args = parser.parse_args(cli)
    return args


def main(cli=None):
    args = parse_args(cli)
    if args.outfile_name is None:
        outfile = os.path.basename(args.snap_name) + '.mat'
    else:
        outfile = args.outfile_name
    convert_to_mat_all_info(args.snap_name, args.density_threshold, outfile)

if __name__ == '__main__':
    main()
    # snap_name = "snaps/snapshot_0149"
    # density_threshold = 6e-6
    # convert_snapshot_to_mat(snap_name, density_threshold, 'snapshot_0149_id.mat')