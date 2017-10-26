import numpy as np
import pynbody
import collections
import tqdm

def compute_cog(snapshots, directory='', save_cache=True, cache_file="cog.npz"):
    if not isinstance(snapshots, collections.Iterable):
        snapshots = [snapshots]
#     print("Processing {} files".format(len(snapshots)))

    cog = np.zeros((3, len(snapshots)), dtype=float)
    times = np.zeros(len(snapshots), dtype=float)
    for i, snap in enumerate(snapshots):
        # duck-typing in case the snapshots are strings
        if not isinstance(snap, pynbody.gadget.GadgetSnap):
            sim = pynbody.load(os.path.join(directory, snap))
        else:
            sim = snap

        print("{:03d} Analysing {} (time {:.2f} Gyr)".format(i, sim.filename, sim.properties['time'].in_units('Gyr')))
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
