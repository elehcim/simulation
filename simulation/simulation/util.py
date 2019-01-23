import os
import glob
import contextlib
import numpy as np

@contextlib.contextmanager
def np_printoptions(*args, **kwargs):
    original = np.get_printoptions()
    np.set_printoptions(*args, **kwargs)
    try:
        yield
    finally:
        np.set_printoptions(**original)


def snapshot_list(dirname, stem="snapshot_", fillwidth=4, include_dir=False):
    """

    Parameters
    ----------
    dirname : str
        Directory with snapshots.
    stem : str
        How snapshot file names start.
    fillwidth : int
        number of the snapshot padded with `fillwidth` zeros.
    include_dir :
        Include the full directory in the output and not the snapshot name only.
    Returns
    -------
    A list of the snapshot files
    """
    if not os.path.isdir(dirname):
        raise IOError("{} is not a directory".format(dirname))
    if include_dir:
        filelist = glob.glob(os.path.join(dirname, stem) + "*")
    else:
        filelist = list(map(os.path.basename, glob.glob(os.path.join(dirname, stem) + "*")))
    filelist.sort()
    return filelist


def first_last_snap(dirname, stem="snapshot_", fillwidth=4):
    if not os.path.isdir(dirname):
        raise IOError("{} is not a directory".format(dirname))
    filelist = snapshot_list(dirname=dirname, stem=stem, include_dir=False)
    first = int(filelist[0][len(stem):])
    last = int(filelist[-1][len(stem):])

    return first, last


def get_snapshot_data(simulation, snap=None):
    import chyplot
    dr = chyplot.CDataGadget()
    fdir = os.path.expanduser(simulation)
    print("Using snapshots in {}".format(fdir))

    dr.setPrefix( fdir )

    if snap is None:
        first_snap, last_snap = first_last_snap(fdir)
        print("Found snapshots [{}: {}]".format(first_snap, last_snap))
        my_snap = last_snap
    else:
        my_snap = snap

    print("Using snapshot {}".format(my_snap))
    dr.set_file(my_snap)

    # dr.checkFilesPresent() # set the first and last dump
    # print dr.lastDump()
    # dr.set_file( dr.lastDump() )

    data = dr.readFile()
    return data