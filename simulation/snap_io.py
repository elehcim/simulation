import os
import glob
import pynbody
import numpy as np
from .util import setup_logger

MORIA_PATH = '/home/michele/sim/MoRIA/'
KICKED_PATH = '/home/michele/sim/MySimulations/Moria8Gyr_tidal'

logger = setup_logger('snap_io', logger_level='INFO')


def make_snaps_path(sim_number, kicked):
    if kicked:
        return os.path.join(KICKED_PATH, 'sim{}'.format(sim_number), 'out')
    else:
        return os.path.join(MORIA_PATH, 'sim{}'.format(sim_number))


def snapshot_file_list(dirname, stem="snapshot_", fillwidth=4, include_dir=False):
    """Return a list of the path to all the snapshots in the simulation folder"""
    if not os.path.isdir(dirname):
         raise IOError("{} is not a directory".format(dirname))
    if include_dir:
         filelist = glob.glob(os.path.join(dirname, stem) + "????")
    else:
         filelist = list(map(os.path.basename, glob.glob(os.path.join(dirname, stem) + "????")))
    filelist.sort()
    return filelist


def load_sim(snap_dir, snap_indexes=None):
    """Return a list of pynbody.SimSnap contained in the directory `snap_dir`"""
    snap_name_list = snapshot_file_list(os.path.expanduser(snap_dir), include_dir=True)
    logger.info("Found {} snapshots".format(len(snap_name_list)))
    if snap_indexes is not None:
        if isinstance(snap_indexes, (list, tuple)):
            snap_name_list = np.array(snap_name_list)[snap_indexes]
        else:
            snap_name_list = snap_name_list[snap_indexes]
        logger.info("Taking {} snapshots ({})".format(len(snap_name_list), snap_indexes))
    snap_list = list(pynbody.load(snap) for snap in snap_name_list)
    return snap_list


def load_snap(snap_dir, snap_number):
    """Return pynbody.SimSnap
    if `snap_number is negative use it as a list index"""
    snaplist = snapshot_file_list(os.path.expanduser(snap_dir), include_dir=True)
    if snap_number < 0:
        return pynbody.load(snaplist[snap_number])
    for snap in snaplist:
        if "{:04d}".format(snap_number) in os.path.basename(snap):
            break
    else:
        raise RuntimeError("No snap {} in simulation folder {}".format(snap_number, snap_dir))
    return pynbody.load(snap)


def load_moria(sim_number, snap_number=None, path=None, snap_indexes=None):
    if path is None:
        sim_dir = make_snaps_path(sim_number, kicked=False)
    else:
        sim_dir = os.path.join(MORIA_PATH, 'sim{}'.format(sim_number))

    if snap_number is None:
        return load_sim(sim_dir, snap_indexes=snap_indexes)
    else:
        return load_snap(sim_dir, snap_number)


def load_kicked(sim_number, snap_number=None, path=None, snap_indexes=None):
    if path is None:
        sim_dir = make_snaps_path(sim_number, kicked=True)
    else:
        sim_dir = os.path.join(KICKED_PATH, 'sim{}'.format(sim_number), 'out')

    if snap_number is None:
        return load_sim(sim_dir, snap_indexes=snap_indexes)
    else:
        return load_snap(sim_dir, snap_number)


def load_moria_sim_and_kicked(sim_number, moria_path=MORIA_PATH, kicked_path=KICKED_PATH):
    return load_moria(sim_number, path=moria_path), load_kicked(sim_number, path=kicked_path)


def load_first_last_snap(snap_dir):
    snaplist = snapshot_file_list(os.path.expanduser(snap_dir), include_dir=True)
    return pynbody.load(snaplist[0]), pynbody.load(snaplist[-1])