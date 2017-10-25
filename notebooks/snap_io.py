import os
import glob
import pynbody

moria_path = '/home/michele/sim/MoRIA/'
kicked_path = '/home/michele/sim/MySimulations/Moria8Gyr_tidal'

def snapshot_list(dirname, stem="snapshot_", fillwidth=4, include_dir=False):
    if not os.path.isdir(dirname):
         raise IOError("{} is not a directory".format(dirname))
    if include_dir:
         filelist = glob.glob(os.path.join(dirname, stem) + "*")
    else:
         filelist = list(map(os.path.basename, glob.glob(os.path.join(dirname, stem) + "*")))
    filelist.sort()
    return filelist

def load_sim(path_to_snapshots):
    snaplist = snapshot_list(path_to_snapshots, include_dir=True)
    simlist = [pynbody.load(snap) for snap in snaplist]
    return simlist

def load_moria_sim_and_kicked(sim_number, kicked=True, moria_path=moria_path, kicked_path=kicked_path):
    snaps_path = os.path.join(moria_path, 'sim{}'.format(sim_number))
    simlist = load_sim(snaps_path)
    
    if kicked:
        ksnaps_path = os.path.join(kicked_path, 'sim{}'.format(sim_number), 'out')
        ksimlist = load_sim(ksnaps_path)
        return simlist, ksimlist
    else:
        return simlist