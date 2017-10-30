import os
import glob
import pynbody

MORIA_PATH = '/home/michele/sim/MoRIA/'
KICKED_PATH = '/home/michele/sim/MySimulations/Moria8Gyr_tidal'

def snapshot_list(dirname, stem="snapshot_", fillwidth=4, include_dir=False):
    """Return a list of the path to all the snapshots in the simulation folder"""
    if not os.path.isdir(dirname):
         raise IOError("{} is not a directory".format(dirname))
    if include_dir:
         filelist = glob.glob(os.path.join(dirname, stem) + "*")
    else:
         filelist = list(map(os.path.basename, glob.glob(os.path.join(dirname, stem) + "*")))
    filelist.sort()
    return filelist

def load_sim(sim_dir):
    """Return a tuple of pynbody.SimSnap contained in the directory `sim_dir`"""
    snaplist = snapshot_list(sim_dir, include_dir=True)
    simlist = tuple(pynbody.load(snap) for snap in snaplist)
    return simlist

def load_snap(sim_dir, snap_number):
    """Return pynbody.SimSnap
    if `snap_number is negative use it as a list index"""
    snaplist = snapshot_list(sim_dir, include_dir=True)
    if snap_number < 0:
        return pynbody.load(snaplist[snap_number])
    for snap in snaplist:
        if "{:04d}".format(snap_number) in os.path.basename(snap):
            break
    else:
        raise RuntimeError("No snap {} in simulation folder {}".format(snap_number, sim_dir))
    return pynbody.load(snap)

def load_moria(sim_number, snap_number=None, path=MORIA_PATH):
    sim_dir = os.path.join(path, 'sim{}'.format(sim_number))
    if snap_number is None:
        return load_sim(sim_dir)
    else:
        return load_snap(sim_dir, snap_number)

def load_kicked(sim_number, snap_number=None, path=KICKED_PATH):
    sim_dir = os.path.join(path, 'sim{}'.format(sim_number), 'out')
    if snap_number is None:
        return load_sim(sim_dir)
    else:
        return load_snap(sim_dir, snap_number)

def load_moria_sim_and_kicked(sim_number, moria_path=MORIA_PATH, kicked_path=KICKED_PATH):
    return load_moria(sim_number, path=moria_path), load_kicked(sim_number, path=kicked_path)