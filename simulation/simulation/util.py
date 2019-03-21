import os
import glob
import contextlib
import numpy as np
import logging
import astropy.units as u

def setup_logger(logger_name=None, logger_level='DEBUG'):

    # formatter = logging.Formatter('%(asctime)s %(levelname)s %(name)s - %(message)s')
    # log_file_handler = logging.FileHandler(os.path.join(config_logger['log_dir'], config_logger['log_file']))
    # log_file_handler.setFormatter(formatter)

    stream_formatter = logging.Formatter('%(asctime)s [%(levelname)-5.5s] - %(message)s')
    log_stream_handler = logging.StreamHandler()
    log_stream_handler.setFormatter(stream_formatter)

    logger = logging.getLogger(logger_name)
    # logger.addHandler(log_file_handler)
    logger.addHandler(log_stream_handler)
    logger.setLevel(logging.getLevelName(logger_level))

    return logger


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


def to_astropy_quantity(simarr, units=None):
    return u.Quantity(simarr.view(type=np.ndarray), unit=units if units is not None else str(simarr.units))


def get_sim_name(sim_path):
    """ Get the name of the simulation

    Example
    -------
    >>> get_sim_name("/home/michele/sim/MySimulations/ng/mb.62002_p200_a800_r600/out")
    mb.62002_p200_a800_r600
    """
    if not(os.path.normpath(sim_path).endswith('out')) and 'MoRIA' in sim_path:
        sim_name = os.path.basename(os.path.normpath(sim_path))
    else:
        sim_name = os.path.basename(os.path.dirname(os.path.normpath(sim_path)))
    return sim_name


def get_sim_traj(sim_name):
    """
    Get sim name and trajectory (+ optional suffix) from the name of the simulation

    Example
    -------
    >>> get_sim_traj(mb.62002_p200_a800_r600_no_gas)
    (mb.62002, p200_a800_r600_no_gas)

    """
    s, t = sim_name.split('002_')
    return s+'002', t


if __name__ == '__main__':
    # mega for-loop
    print(good_sims)

    for sim_name in map(os.path.basename, good_sims):
        SIM, TRAJ = get_sim_traj(sim_name)
        sim_path = os.path.join(SIMPATH, "{}_{}".format(SIM, TRAJ), "out")