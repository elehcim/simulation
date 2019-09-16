import pynbody
import numpy as np
import os
from scipy.io import savemat



def convert_snapshot_to_mat(snap_name, density_threshold, outfile_name):
    """Take a snapshot, cut on the density and provide a mat file with particle positions (x,y,z) and IDs"""

    snap = pynbody.load(snap_name)
    gas = snap.gas
    filt = pynbody.filt.HighPass('rho', density_threshold)

    pos = gas[filt]['pos'].view(np.ndarray)
    ids = gas[filt]['iord'].view(np.ndarray)
    rho = gas[filt]['rho'].view(np.ndarray)

    data = dict(x=pos[:,0], y=pos[:,1], z=pos[:,2], id=ids, rho=rho)
    print('output file', outfile_name)
    savemat(outfile_name, data)
    return pos

if __name__ == '__main__':
    snap_name = "snaps/snapshot_0149"
    density_threshold = 6e-6
    convert_snapshot_to_mat(snap_name, density_threshold, 'snapshot_0149_id.mat')