import os
import glob
import contextlib
import numpy as np
import logging
import json
import astropy.units as u
from astropy.table import Table

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
    return u.Quantity(simarr.view(type=np.ndarray), unit=units if units is not None else str(simarr.units), dtype=simarr.dtype)


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


def contains_duplicates(X):
    return len(np.unique(X)) != len(X)


def make_df_monotonic_again(df, col='t'):
    """Return a copy of the df where the non monotonic parts are cut out"""
    diff = df.t.diff()
    restart_points = df.t[diff <0]
    before_restart_points = df.t.loc[restart_points.index-1]
    end_idx = [df.query('index > {} and t > {}'.format(idx, t)).index[0] for idx, t in zip(before_restart_points.index, before_restart_points.values)]
    na = df.copy()
    for a, b in zip(before_restart_points.index, end_idx):
        na = na.drop(df.index[slice(a,b)])
    return na


def get_no_gti_intervals(info):
    diff = info.step.diff()
    restart_points = info.step[diff < 0]
    idx_restart = [idx - ((info.step.loc[idx - 1]) - v + 1) for idx, v in zip(restart_points.index, restart_points.values)]
    return idx_restart, restart_points


def prune_no_gti(df, idx_restart, restart_points):
    na = df.copy()
    for a, b in zip(idx_restart, restart_points.index):
        na = na.drop(df.index[slice(a,b)])
    return na


def make_info_monotonic_again(info):
    diff = info.step.diff()
    restart_points = info.step[diff < 0]
    idx_restart = [idx - ((info.step.loc[idx - 1]) - v + 1) for idx, v in zip(restart_points.index, restart_points.values)]

    na = info.copy()
    for a, b in zip(idx_restart, restart_points.index):
        na = na.drop(info.index[slice(a,b)])
    return na


def make_df_monotonic_again_using_info(df, info):
    # FIXME not tested
    assert len(df) == len(info), "I can use this function only if df is the same length as info"
    idx_restart, restart_points = get_no_gti_intervals(info)
    new_df = prune_no_gti(df, idx_restart, restart_points)
    return new_df


def get_omega_mb(sim_name, omega_dir='~/sim/analysis/ng_ana/data/quat'):
    """Return a numpy array reading the omega_mb table in `quat_dir`"""
    logger = setup_logger('get_omega_mb', logger_level='INFO')
    if os.path.isdir(os.path.expanduser(omega_dir)):
        omega_dir = os.path.expanduser(omega_dir)
        omega_file = os.path.join(omega_dir, sim_name+'_quat.fits')
    else:
        omega_file = None

    if os.path.isfile(omega_file):
        logger.info('Reading omega_mb table: {}'.format(omega_file))
        tbl = Table.read(omega_file)
        omega_mb_arr = np.array([tbl["omega_mb_x"], tbl["omega_mb_y"], tbl["omega_mb_z"]]).T
    else:
        logger.warning('Cannot find omega_mb table...')
        omega_mb_arr = None
    return omega_mb_arr

def get_quat(sim_name, quat_dir='~/sim/analysis/ng_ana/data/quat'):
    """Return a numpy array reading the table in `quat_dir`"""
    logger = setup_logger('get_quat', logger_level='INFO')
    if os.path.isdir(os.path.expanduser(quat_dir)):
        quat_dir = os.path.expanduser(quat_dir)
        quat_file = os.path.join(quat_dir, sim_name+'_quat.fits')
    else:
        quat_file = None

    if os.path.isfile(quat_file):
        logger.info('Reading quaternion table: {}'.format(quat_file))
        tbl = Table.read(quat_file)
        quat_arr = np.array([tbl["q_w"], tbl["q_x"], tbl["q_y"], tbl["q_z"]]).T
    else:
        logger.warning('Cannot find quaternion table...')
        quat_arr = None
    return quat_arr


def get_pivot(sim_name,
              pivot_file='~/sim/MySimulations/ng/pivot.json',
              raise_if_cannot_derotate=True):
    logger = setup_logger('get_pivot', logger_level='INFO')
    try:
        with open(os.path.expanduser(pivot_file), 'r') as f:
            d = json.load(f)['pivot']
        pivot = np.array(d[sim_name].split(), dtype=np.float64)
    except Exception as e:
        if raise_if_cannot_derotate:
            raise e
        else:
            print(e)
            logger.warning('Cannot find pivot table...')
            pivot = None
    return pivot


def get_quat_pivot(sim_name,
                   quat_dir='~/sim/analysis/ng_ana/data/quat',
                   pivot_file='~/sim/MySimulations/ng/pivot.json',
                   raise_if_cannot_derotate=True):

    quat_arr = get_quat(sim_name, quat_dir)
    pivot = get_pivot(sim_name, pivot_file, raise_if_cannot_derotate)
    return quat_arr, pivot


def get_quat_omega_pivot(sim_name,
                   quat_dir='~/sim/analysis/ng_ana/data/quat',
                   pivot_file='~/sim/MySimulations/ng/pivot.json',
                   raise_if_cannot_derotate=True):

    quat_arr = get_quat(sim_name, quat_dir)
    omega_mb_arr = get_omega_mb(sim_name, quat_dir)
    pivot = get_pivot(sim_name, pivot_file, raise_if_cannot_derotate)
    return quat_arr, omega_mb_arr, pivot


if __name__ == '__main__':
    # mega for-loop
    print(good_sims)

    for sim_name in map(os.path.basename, good_sims):
        SIM, TRAJ = get_sim_traj(sim_name)
        sim_path = os.path.join(SIMPATH, "{}_{}".format(SIM, TRAJ), "out")