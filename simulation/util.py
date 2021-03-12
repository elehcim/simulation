import os
import glob
import contextlib
import numpy as np
import pandas as pd
import logging
import json
import astropy.units as u
from astropy.table import Table

DATA_DIR = os.getenv('SIM_DATA_DIR', default='/home/michele/sim/analysis/ng_ana/data')
SIMS_DIR = os.getenv('SIM_SIM_DIR', default='/home/michele/sim/MySimulations/ok_new_adhoc_or_not_affected')

PIVOT_FILE = os.path.join(os.path.dirname(__file__), 'pivot.json')
DEROTATION_DIR = os.path.join(DATA_DIR, 'quat')

loggers = {}


def setup_logger(logger_name=None, logger_level='DEBUG'):
    if loggers.get(logger_name):
        logger = loggers.get(logger_name)
    else:


    # formatter = logging.Formatter('%(asctime)s %(levelname)s %(name)s - %(message)s')
    # log_file_handler = logging.FileHandler(os.path.join(config_logger['log_dir'], config_logger['log_file']))
    # log_file_handler.setFormatter(formatter)

        stream_formatter = logging.Formatter('%(asctime)s (%(name)s) [%(levelname)-5.5s] - %(message)s')
        log_stream_handler = logging.StreamHandler()
        log_stream_handler.setFormatter(stream_formatter)

        logger = logging.getLogger(logger_name)
        # logger.addHandler(log_file_handler)
        logger.addHandler(log_stream_handler)
        logger.setLevel(logging.getLevelName(logger_level))
        loggers[logger_name] = logger

    return logger

logger = setup_logger('derotation_info', logger_level='INFO')


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
    """
    Get intervals which can then be eliminated with something like:

    for a, b in zip(idx_restart, restart_points.index):
        new_df = new_df.drop(df.index[slice(a,b)])

    Return:
    ------
    idx_restart: list
        left index of the intervals to be removed.
    restart_points: pd.Series
        the index is the right index of the intervals to be removed.
        are the value of the index in correspondance of the first value of the restarted run

    Example:
    -------
    >>> test_info = pd.DataFrame({'step': [0, 1, 2, 1, 2, 3, 4, 5, 3, 4, 5, 6, 7]})
    >>> idx_restart, restart_point = get_no_gti_intervals(test_info)
    >>> idx_restart
        [1, 5]
    >>> restart_point.index
        [3, 8]

    # To have a monotonic I should remove from index 1 to index 3 and then from 5 to 8 (keeping the right one)

      idx:          0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12
      step:         0, 1, 2, 1, 2, 3, 4, 5, 3, 4, 5, 6, 7
      to be removed:   |--|<       |-----|<

      result:       0,       1, 2,          3, 4, 5, 6, 7

    # After calling also `prune_no_gti`
      idx:          0,       3, 4,          8, 9,10,11,12
      step:         0,       1, 2,          3, 4, 5, 6, 7
    """
    diff = info.step.diff()
    restart_points = info.step[diff < 0]
    idx_restart = [idx - ((info.step.loc[idx - 1]) - v + 1) for idx, v in zip(restart_points.index, restart_points.values)]
    return idx_restart, restart_points


def prune_no_gti(df, idx_restart, restart_points):
    na = df.copy()
    for a, b in zip(idx_restart, restart_points.index):
        na = na.drop(df.index[slice(a,b)])
    return na.reset_index(drop=True)


def make_info_monotonic_again(info):
    diff = info.step.diff()
    restart_points = info.step[diff < 0]
    idx_restart = [idx - ((info.step.loc[idx - 1]) - v + 1) for idx, v in zip(restart_points.index, restart_points.values)]

    na = info.copy()
    for a, b in zip(idx_restart, restart_points.index):
        na = na.drop(info.index[slice(a,b)])
    return na


def make_df_monotonic_again_using_reference_df(df, ref):
    # FIXME not tested
    assert len(df) == len(ref), "I can use this function only if df is the same length as reference df (df:{}, ref:{})".format(len(df), len(ref))
    idx_restart, restart_points = get_no_gti_intervals(ref)
    new_df = prune_no_gti(df, idx_restart, restart_points)
    return new_df


def get_omega_mb(sim_name, omega_dir=DEROTATION_DIR):
    """Return a numpy array reading the omega_mb table in `quat_dir`"""
    if os.path.isdir(os.path.expanduser(omega_dir)):
        omega_dir = os.path.expanduser(omega_dir)
        omega_file = os.path.join(omega_dir, sim_name+'_quat.fits')
    else:
        omega_file = None

    if os.path.isfile(omega_file):
        logger.debug('Reading omega_mb table: {}'.format(omega_file))
        tbl = Table.read(omega_file)
        omega_mb_arr = np.array([tbl["omega_mb_x"], tbl["omega_mb_y"], tbl["omega_mb_z"]]).T
    else:
        logger.warning('Cannot find omega_mb table...')
        omega_mb_arr = None
    return omega_mb_arr


def get_quat_file(sim_name, quat_dir=DEROTATION_DIR):
    if os.path.isdir(os.path.expanduser(quat_dir)):
        quat_dir = os.path.expanduser(quat_dir)
        quat_file = os.path.join(quat_dir, sim_name+'_quat.fits')
    else:
        quat_file = None
    return quat_file


def get_quat(sim_name, quat_dir=DEROTATION_DIR):
    """Return a numpy array reading the table in `quat_dir`"""
    quat_file = get_quat_file(sim_name, quat_dir)
    if os.path.isfile(quat_file):
        logger.debug('Reading quaternion table: {}'.format(quat_file))
        tbl = Table.read(quat_file)
        quat_arr = np.array([tbl["q_w"], tbl["q_x"], tbl["q_y"], tbl["q_z"]]).T
    else:
        logger.warning('Cannot find quaternion table...')
        quat_arr = None
    return quat_arr


def get_initial_rotation(sim_name):
    """Get initial rotation due to the choice of the axis in the moving box technique,
    Depends on the pericenter distance only"""
    quat_file = get_quat_file(sim_name)
    tbl = Table.read(quat_file)
    quat_vp0 = np.quaternion(tbl.meta['QUAT_P0W'], tbl.meta["QUAT_P0X"], tbl.meta["QUAT_P0Y"], tbl.meta["QUAT_P0Z"])
    return quat_vp0


def get_pivot(sim_name,
              pivot_file=PIVOT_FILE,
              raise_if_cannot_derotate=True):
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
                   quat_dir=DEROTATION_DIR,
                   pivot_file=PIVOT_FILE,
                   raise_if_cannot_derotate=True):

    quat_arr = get_quat(sim_name, quat_dir)
    pivot = get_pivot(sim_name, pivot_file, raise_if_cannot_derotate)
    return quat_arr, pivot


def get_quat_omega_pivot(sim_name,
                   quat_dir=DEROTATION_DIR,
                   pivot_file=PIVOT_FILE,
                   raise_if_cannot_derotate=True):

    quat_arr = get_quat(sim_name, quat_dir)
    omega_mb_arr = get_omega_mb(sim_name, quat_dir)
    if quat_arr is None or omega_mb_arr is None:
        if raise_if_cannot_derotate:
            raise RuntimeError('Cannot derotate')
        else:
            logger.warning("I don't have all the ingredients to derotate...")
    pivot = get_pivot(sim_name, pivot_file, raise_if_cannot_derotate)

    return quat_arr, omega_mb_arr, pivot



def make_lowess(series, **kwargs):
    """Adapted from here: https://www.allendowney.com/blog/2019/04/01/local-regression-in-python/"""
    from statsmodels.nonparametric.smoothers_lowess import lowess
    endog = series.values
    exog = series.index.values

    smooth = lowess(endog, exog, **kwargs)
    index, data = np.transpose(smooth)

    return pd.Series(data, index=index)


def parse_simname(name):
    i = name.find('_p')
    j = name.find('_a')
    k = name.find('_r')
    p  = int(name[i+2:j])
    a  = int(name[j+2:k])
    r  = int(name[k+2:k+5])  # Assuming 3 digits
    return p, a, r


def get_aspect(extent, shape):
    """When importing an image with imread, this function computes the aspect parameter.

    Example:
    --------
    img = plt.imread("Mr_gr_-19_-8_0.0_1.0.png")
    extent = (-19, -8, 0, 1)
    ax.imshow(img, extent=extent, aspect=get_aspect(extent, img.shape)

    """
    W = np.abs(extent[1]-extent[0])
    H = np.abs(extent[3]-extent[2])
    h = shape[0]
    w = shape[1]
    return (h/H) / (w/W)

def savefig(fig, file_stem, ext, dpi=300, tight=True, **kwargs):
    file_name = file_stem + ext
    print(f'Saving {file_name}...')
    if tight:
        kwargs.update(bbox_inches='tight')
    print(kwargs)
    if ext == '.png':
        fig.savefig(file_name, dpi=dpi, **kwargs)
        out = f"{file_stem}-crop.png"
        os.system(f'convert -trim {file_name} {out}')
    elif ext == '.pdf':
        fig.savefig(file_name, dpi=dpi, **kwargs)
        os.system(f'pdfcrop {file_name}')

if __name__ == '__main__':
    # mega for-loop
    print(good_sims)

    for sim_name in map(os.path.basename, good_sims):
        SIM, TRAJ = get_sim_traj(sim_name)
        sim_path = os.path.join(SIMPATH, "{}_{}".format(SIM, TRAJ), "out")