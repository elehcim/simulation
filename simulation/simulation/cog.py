import os
import numpy as np
import pynbody
import collections
import tqdm

def compute_cog(snapshots, directory='', save_cache=False, cache_file="cog.npz", verbose=True, family=None):
    """
    Compute the center of gravity of a simulation

    Args:
        snapshots: Iterable of snapshots simulation. Can be either an pynbody.snapshot.SimSnap or a list of string (filenames).
            In the latter case the `directory` argument is used as the containing folder.
        save_cache: Save the resultsin a npz file
        cache_file: Name of the cache npz file 
        verbose: Print which file is currently read
        family: The subfamily attribute to take into account, if `None` consider the whole snapshot.

    Returns:
        Tuple of numpy arrays (t, cog) where:
            `t` is the array of the times of the snapshots in Gyr
            `cog` is a 3 column array with the coordinates of the center of gravity positions for each snapshot
    """
    if not isinstance(snapshots, collections.Iterable):
        snapshots = [snapshots]

    cog = np.zeros((3, len(snapshots)), dtype=float)
    times = np.zeros(len(snapshots), dtype=float)
    for i, snap in enumerate(snapshots):
        # duck-typing in case the snapshots are strings
        if not isinstance(snap, (pynbody.gadget.GadgetSnap, pynbody.snapshot.FamilySubSnap)):
            sim = pynbody.load(os.path.join(directory, snap))
        else:
            sim = snap

        if family is not None:
            sim = sim.__getattr__(family)

        if verbose:
            print("{:03d} Analysing {} (time {:.4f} Gyr)".format(i, sim.filename, sim.properties['time'].in_units('Gyr')))
        times[i] = sim.properties['time'].in_units('Gyr')
        mass = sim['mass']
        pos = sim['pos']

        tot_mass = mass.sum()
        # print "tot mass = {} Msol".format(tot_mass)
        # cog[:,i] = np.array([(pos[:,0] * mass).sum(), (pos[:,1] * mass).sum(), (pos[:,2] * mass).sum()])/tot_mass
        cog[:,i] = np.sum(mass * pos.transpose(), axis=1) / tot_mass
    if save_cache:
        np.savez(cache_file, times=times, cog=cog)
    return times, cog
