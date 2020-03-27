import os
import simulation
import pandas as pd
import matplotlib.pylab as plt
import numpy as np
import tqdm
import glob
import pickle
import warnings
from astropy.table import Table
from simulation.util import make_lowess, get_sim_name, setup_logger, DATA_DIR, SIMS_DIR
from simulation.derived import feh, mgfe, gas_metals, neutral_fraction
from simulation.vlos_profiles import get_max_vlos
from simulation.magnitude_transformation import get_sdss_u, get_sdss_g, get_sdss_r, get_sdss_i, get_sdss_z

logger = setup_logger('simdata', logger_level='INFO')
# TABLE_LIST_DIR = '/home/michele/sim/analysis/ng_ana/data/tables/general/'

SIM_NAME_DICT = {
 '62p50': 'mb.62002_p50_a800_r600',
 '62p100': 'mb.62002_p100_a800_r600',
 '62p150': 'mb.62002_p150_a800_r600',
 '62p200': 'mb.62002_p200_a800_r600',
 '62p300': 'mb.62002_p300_a800_r600',
 '71p50': 'mb.71002_p50_a800_r600',
 '71p100': 'mb.71002_p100_a800_r600',
 '71p150': 'mb.71002_p150_a800_r600',
 '71p200': 'mb.71002_p200_a800_r600',
 '71p300': 'mb.71002_p300_a800_r600',
 '68p50': 'mb.68002_p50_a800_r600',
 '68p100': 'mb.68002_p100_a800_r600',
 '68p150': 'mb.68002_p150_a800_r600',
 '68p200': 'mb.68002_p200_a800_r600',
 '68p300': 'mb.68002_p300_a800_r600',
 '69p50': 'mb.69002_p50_a800_r600_new_adhoc',
 '69p100': 'mb.69002_p100_a800_r600_new_adhoc',
 '69p150': 'mb.69002_p150_a800_r600_new_adhoc',
 '69p200': 'mb.69002_p200_a800_r600_new_adhoc',
 '69p300': 'mb.69002_p300_a800_r600_new_adhoc',
 '41p50': 'mb.41002_p50_a800_r600_t9.56',
 '41p100': 'mb.41002_p100_a800_r600_t9.56',
 '41p150': 'mb.41002_p150_a800_r600_t9.56',
 '41p200': 'mb.41002_p200_a800_r600_t9.56',
 '41p300': 'mb.41002_p300_a800_r600_t9.56',
}

NO_GAS_DICT = {
 '62p50': 'mb.62002_p50_a800_r600_no_gas',
 '62p100': 'mb.62002_p100_a800_r600_no_gas',
 '62p150': 'mb.62002_p150_a800_r600_no_gas',
 '62p200': 'mb.62002_p200_a800_r600_no_gas',
 '62p300': 'mb.62002_p300_a800_r600_no_gas',
 '69p50': 'mb.69002_p50_a800_r600_no_gas',
 '69p100': 'mb.69002_p100_a800_r600_no_gas',
 '69p150': 'mb.69002_p150_a800_r600_no_gas',
 '69p200': 'mb.69002_p200_a800_r600_no_gas',
 '69p300': 'mb.69002_p300_a800_r600_no_gas',
 '41p50': 'mb.41002_p50_a800_r600_t9.56_no_gas',
 '41p100': 'mb.41002_p100_a800_r600_t9.56_no_gas',
 '41p150': 'mb.41002_p150_a800_r600_t9.56_no_gas',
 '41p200': 'mb.41002_p200_a800_r600_t9.56_no_gas',
 '41p300': 'mb.41002_p300_a800_r600_t9.56_no_gas',
}

SPECIAL_DICT = {
 '68p100r': 'mb.68002_p100_a800_r600_retrograde',
 '69p200i': 'mb.69002_p200_a800_r600_inertial',

}


def get_name_peri(sim_name):
    mb, peri = sim_name.split('_')[:2]
    return mb[3:], int(peri[1:])


def get_radial_period_pickle(sim_name, data_dir=DATA_DIR):
    """Return radial period in Gyr"""
    mb, peri = sim_name.split('_')[:2]
    filename = os.path.join(data_dir, "radial_period", "radial_period_{}.pickle".format(mb))
    # print(filename)
    d = pickle.load(open(filename, 'rb'))
    # print(d)
    rperiod = d[int(peri[1:])]
    return rperiod


def get_radial_period(sim_name, which='measured', data_dir=DATA_DIR):
    """Return radial period in Gyr

    Parameters
    ----------
    which : str
        one of ['measured', 'computed'] the difference being that the latter has
        been computing using a potential for the cluster and the information of
        the orbit pericenter and apocenter.
    """

    name = os.path.join(data_dir, f"radial_period/{sim_name}_radial_period.fits")
    logger.debug(f"Getting T_r table: {name}")
    return float(Table.read(name)[which])


def compute_t_period(sim_name, df=None):
    """Basically return `tables` with some additional useful orbit related columns.
    You can specify the version (orbit_sideon or not) of the table through parameter `df`.
    But this does not change the result of the added columns because they are independent from the point of view.
    """

    if df is None:
        df = get_tables(sim_name, orbit_sideon=True).to_pandas()

    # Some have it, some not
    if 'r' not in df.keys():
        df['r'] = np.sqrt(df['x']**2 + df['y']**2)

    radial_period = get_radial_period(sim_name)
    zero_crossings = np.where(np.diff(np.signbit(df['r'].diff())))[0]
    idx_peri = zero_crossings[1]  # first element is always 0, so the second is the actual pericenter
    first_pericenter_time = df['t'][idx_peri]
    # print(zero_crossings)
    # print(first_pericenter_time)

    # Put 0 of the scale on first pericenter
    df['t_period'] = (df['t']-first_pericenter_time)/radial_period
    df['orbital_phase'] = pd.cut(df['t_period'], [-0.25, 0.25, 0.75, 1.25, 1.75, 2.25], labels=False)
    df['offset_orbital_phase'] = pd.cut(df['t_period'], np.array([-0.25, 0.25, 0.75, 1.25, 1.75, 2.25])-0.25, labels=False)
    return df


def get_t_period(sim_name):
    return compute_t_period(sim_name)['t_period']


def get_vrot_max(tbl):
    max_vrot = list()
    for prof in tbl['v_circ']:
        max_vrot.append(np.nanmax(prof))
    return np.array(max_vrot)


def get_center(sim_name, data_dir=DATA_DIR):
    tbl = Table.read(os.path.join(data_dir, 'center', sim_name + '_cen.fits'))
    cen = np.array([tbl['cx'], tbl['cy'], tbl['cz']]).T
    return cen


def get_df(sim_name, window_size=20, std=30, cut=None, data_dir=DATA_DIR):
    # TODO improve API with explicit `orbit_sideon` parameter
    logger.info("Merging dataframes")
    df = Table.read(os.path.join(data_dir, 'tables/{}.fits'.format(sim_name))).to_pandas()

    # Remove wrong sigmas
    unwanted_columns = ['sigma_star', 'sigma_gas']
    for u in unwanted_columns:
        if u in df.columns: del df[u]

    is_sideon = is_orbit_sideon(sim_name)
    name_no_orientation = sim_name.split('_orbit_sideon')[0]

    phot_tbl = get_phot(name_no_orientation, orbit_sideon=is_sideon, data_dir=data_dir).to_pandas()

    struct_tbl = get_structure(name_no_orientation, orbit_sideon=is_sideon).to_pandas()
    dm_tbl = get_dm(name_no_orientation).to_pandas()
    rt_tbl = get_tidal_radius(name_no_orientation).to_pandas()
    vr_tbl = get_virial(name_no_orientation).to_pandas()
    mag_tbl = get_magnitudes(name_no_orientation).to_pandas()
    cg_tbl = get_cold_gas(name_no_orientation).to_pandas()
    sf_tbl = get_sf(name_no_orientation).to_pandas()
    sig_tbl = get_sigma(name_no_orientation, orbit_sideon=is_sideon).to_pandas()
    lr_tbl = get_lambda_r(name_no_orientation, orbit_sideon=is_sideon)
    am_tbl = get_angmom(name_no_orientation, orbit_sideon=is_sideon)

    # Merge data
    for col in phot_tbl.columns:
        # skip this:
        if col == 'lambda_r': continue
        df[col] = phot_tbl[col]

    for col in dm_tbl.columns:
        df[col] = dm_tbl[col]

    for col in rt_tbl.columns:
        df[col] = rt_tbl[col]

    for col in vr_tbl.columns:
        df[col] = vr_tbl[col]

    for col in sig_tbl.columns:
        df[col] = sig_tbl[col]

    for col in cg_tbl.columns:
        df[col] = cg_tbl[col]

    for col in mag_tbl.columns:
        df['mag_' + col] = mag_tbl[col]

    for col in ('n', 'mu_e', 'avg_mu_e'):
        df[col] = struct_tbl[col]

    for col in ('medians', 'averages', 'length', 'bx', 'by', 'bz'):
        df['sf_' + col] = sf_tbl[col]

    for col in ['lambda_r']:
        df[col] = lr_tbl[col]

    for col in am_tbl.columns:
        # I already have those:
        if col in ('r', 't', 't_period', 'orbital_phase', 'mass_star'):
            continue
        df[col] = am_tbl[col]

    for c in ('rms_err', 'exit_mode', 'numiter', 'r_eff'):
        df['fit_' + c] = struct_tbl[c]

    df['name'], df['pericenter'] = get_name_peri(sim_name)


    avg_columns = ['lambda_r', 'n', 'sfr', 'r_eff', 'r_eff3d',
                   'fit_r_eff', 'ellipticity', 'mu_e', 'avg_mu_e',
                   'sigma_gas', 'sigma_star',
                   'mass_star', 'orientation', 'dm_mass', 'dm_mass_all']

    # Max v_los
    try:
        df['max_v_los'] = get_max_vlos(sim_name, is_sideon)
        avg_columns.append('max_v_los')
    # FIXME ok but this is not what I need
    except Exception as e:
        logger.error(e)
        logger.error("Not computing Max vlos")

    # Do averages

    for col in avg_columns:
        # Do not use spashot with no central stars to compute lowes on sigma.
        if col.startswith('sigma'):
            df[col].replace(0.0, np.nan, inplace=True)

        # LOWESS
        df[col+'_mean'] = make_lowess(df[col], frac=0.1)

        # df[col+'_mean'] = df[col].rolling(window_size, center=True).mean()
        df[col+'_std'] = df[col].rolling(window_size).std()
        if col.startswith('sigma'):
            df[col+'_mean'][df[col+'_mean'] <= 0.0] = np.nan
        # Can use bfill()) or min_periods
        # df[col+'_mean'] = df[col].rolling(window=window_size, win_type='gaussian', center=True, min_periods=1).mean(std=std)
        # print(col)
        # print(df[col].dtype.byteorder)
        # if df[col].dtype.byteorder == '>':  # big-endian, i.e. non native
        #     ser = pd.Series(df[col].to_numpy().byteswap().newbyteorder())
        # else:
        #     ser = df[col]
        # df[col+'_mean'] = ser.rolling(10).mean()

    # Cut if needed, (Moria)
    if cut is not None:
        df = df.iloc[cut:].copy()

    df['r'] = np.sqrt(df.x**2 + df.y**2)

    if 'max_v_los' in df.columns and 'sigma_star' in df.columns:
        df['v_over_sigma'] = df['max_v_los']/df['sigma_star']
        # LOWESS
        df['v_over_sigma'+'_mean'] = make_lowess(df['v_over_sigma'], frac=0.1)
        # STD formula in case of division https://en.wikipedia.org/wiki/Propagation_of_uncertainty#Example_formulae
        df['v_over_sigma_std'] = df['v_over_sigma']*np.sqrt((df['max_v_los_std']/df['max_v_los'])**2 + (df['sigma_star' + '_std']/df['sigma_star'])**2)

    # t_period
    try:
        radial_period = get_radial_period(name_no_orientation, data_dir=data_dir)
        zero_crossings = np.where(np.diff(np.signbit(df['r'].diff())))[0]
        idx_peri = zero_crossings[1]  # first element is always 0, so the second is the actual pericenter
        first_pericenter_time = df.t[idx_peri]
        # print(zero_crossings)
        # print(first_pericenter_time)

        # Put 0 of the scale on first pericenter
        df['t_period'] = (df.t-first_pericenter_time)/radial_period
        t_period_intervals = np.array([-0.25, 0.25, 0.75, 1.25, 1.75, 2.25])

        df['orbital_phase'] = pd.cut(df.t_period, t_period_intervals, labels=False)
        df['offset_orbital_phase'] = pd.cut(df.t_period, t_period_intervals - 0.25, labels=False)
            # labels=(1, 2, 3, 4))
            # labels=('peri_1', 'apo_1', 'peri_2', 'apo_2' ))
        # 0 on start of simulation
        # df['t_period'] = (df.t-df.t.iloc[0])/radial_period

    except Exception as e:
        logger.error(repr(e), exc_info=True)
        logger.error("Not computing Radial period")

    logger.info('{}: ok, data loaded'.format(sim_name))
    return df

def get_last_d(d):
    last_d = dict()
    for k, df in d.items():
        last_row = get_df_last_rows(df)
        last_d[k] = last_row
    return last_d

def get_df_last_rows(df, n_last=5):
    """From a dataframe get the last row among the last n_last if they are available."""
    last_row = pd.DataFrame(data=None, columns=df.columns, index=(0,))
    for c in df.columns:
        lvi = df[c].last_valid_index()
        for i in range(1, n_last + 1):
            if lvi == len(df) - i: # if it's one of the n_last rows it's still ok
                last_row[c] = df[c].iloc[lvi]
                break
        else:
            if not (c.endswith('_mean') or c.endswith('_std') or c.startswith('sf_')):
                logger.debug(f'In {df.name[0][:3]}p{df.pericenter[0]}tf{df.t.iloc[-1]:.2f} {c} we have a different number of elements ({lvi}, {df.t.iloc[lvi]:.2f})')
            last_row[c] = np.nan
    return last_row

def get_first_d(d):
    first_d = dict()
    for k, df in d.items():
        first_d[k] = df.iloc[[0]]  # double brackets return Dataframe, not Series
    return first_d

def last_d2df(last_d):
    last_df = pd.concat(last_d).reset_index(level=1, drop=True)
    # last_df.columns

    drop_cols = ['bbox_xmax', 'bbox_xmin', 'bbox_ymax', 'bbox_ymin', 'area',
                 'xcentroid', 'xmax','xmin', 'ycentroid', 'ymax','ymin',
                 'cxx', 'cyy', 'cxy', 'covar_sigx2', 'covar_sigxy', 'covar_sigy2',
                 'gini', 't', 'perimeter', 'min_value', 'minval_xpos', 'minval_ypos',
                 'semimajor_axis_sigma', 'semiminor_axis_sigma', 'source_sum',
                 'max_value', 'maxval_xpos', 'maxval_ypos','id','equivalent_radius',
                ]

    not_wanted_cols = list()
    for c in last_df.columns:
        if c in drop_cols or c.endswith('_mean') or c.endswith('_std') or c.startswith('sf_'):
            # print(c)
            not_wanted_cols.append(c)
    last_df.drop(columns=not_wanted_cols, inplace=True)
    # last_df.columns

    # last_df['color_sdss_g_r'] = color_sdss_g_r(last_df['mag_v'] - last_df['mag_r'])
    last_df['color_sdss_g_r'] = last_df['mag_sdss_g'] - last_df['mag_sdss_r']

    # Categorize name column
    last_df['name'] = pd.Categorical(last_df.name, ['62002', '71002', '69002', '68002', '41002'], ordered=True)
    # last_df.name

    last_df.sort_values(by=['name', 'pericenter'], inplace=True)
    return last_df

def get_sf_pos(sim_name):
    return pickle.load(open(os.path.join(DATA_DIR, 'sf', sim_name + "_sf_pos.pkl"), 'rb'))

def get_sf(sim_name, data_dir=DATA_DIR):
    name = os.path.join(data_dir, f"sf/{sim_name}_sf.fits")
    logger.debug(f"Getting Star Formation table: {name}")
    return Table.read(name)

def get_mach(sim_name, data_dir=DATA_DIR):
    name = os.path.join(data_dir, f"mach/{sim_name}_mach.fits")
    logger.debug(f"Getting Mach number table: {name}")
    return Table.read(name)


def shorten_name(filename):
    name, peri = os.path.basename(filename).split('_')[:2]
    return name[3:5] + peri


def get_color_from_name(name):
    prop_cycle = plt.rcParams['axes.prop_cycle']
    mpl_colors = prop_cycle.by_key()['color']
    if name.startswith('41'):
        col = mpl_colors[0]
    elif name.startswith('68'):
        col = mpl_colors[1]
    elif name.startswith('69'):
        col = mpl_colors[2]
    elif name.startswith('71'):
        col = mpl_colors[3]
    elif name.startswith('62'):
        col = mpl_colors[4]
    else:
        col = None
    return col

LINESTYLES_DICT = dict((('loosely dotted',        (0, (1, 10))),
     ('dotted',                (0, (1, 1))),
     ('densely dotted',        (0, (1, 1))),

     ('loosely dashed',        (0, (5, 10))),
     ('dashed',                (0, (5, 5))),
     ('densely dashed',        (0, (5, 1))),

     ('loosely dashdotted',    (0, (3, 10, 1, 10))),
     ('dashdotted',            (0, (3, 5, 1, 5))),
     ('densely dashdotted',    (0, (3, 1, 1, 1))),

     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
     ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
     ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))))


def get_styles_from_name(*args, **kwargs):
    warnings.warn("This function is deprecaate, use simulation.simdata.get_styles_from_peri", DeprecationWarning)
    return get_styles_from_peri(*args, **kwargs)

def get_styles_from_peri(name, scatter=False):
    if scatter:
        styles = ['o', '^', "s", 'd', '*']
    else:
        styles = ['-', '--', '-.', ':', LINESTYLES_DICT['densely dashdotdotted']]
    if name.endswith('p50'):
        st = styles[0]
    elif name.endswith('100'):
        st = styles[1]
    elif name.endswith('150'):
        st = styles[2]
    elif name.endswith('200'):
        st = styles[3]
    elif name.endswith('300'):
        st = styles[4]
    else:
        st = None
    return st



def get_styles_from_sim(name, scatter=False):
    if scatter:
        styles = ['o', '^', "s", 'd', '*']
    else:
        styles = ['-', '--', '-.', ':', LINESTYLES_DICT['densely dashdotdotted']]
    if name.startswith('41'):
        st = styles[0]
    elif name.startswith('69'):
        st = styles[1]
    elif name.startswith('68'):
        st = styles[2]
    elif name.startswith('71'):
        st = styles[3]
    elif name.startswith('62'):
        st = styles[4]
    else:
        st = None
    return st

def get_color_styles(d, scatter=False):
    col_d = dict()
    for k in d.keys():
        col_d[k] = get_color_from_name(k)
    st_d = dict()
    for k in d.keys():
        st_d[k] = get_styles_from_peri(k, scatter=scatter)
    return col_d, st_d

# From https://www.tutorialspoint.com/How-to-correctly-sort-a-string-with-a-number-inside-in-Python
import re

def _atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(_):
    warnings.warn("This function is deprecaate, use simulation.simdata.SIM_NAME_DICT", DeprecationWarning)
    return _natural_keys(_)

def _natural_keys(text):
    return [ _atoi(c) for c in re.split(r'(\d+)',text) ]


def is_orbit_sideon(sim_name):
    if "_orbit_sideon" in sim_name:
        return True
    else:
        return False


# Legacy function
# def get_sim_name_list(table_list_dir=TABLE_LIST_DIR):
#     tables_list = glob.glob(os.path.join(table_list_dir, '*.fits'))
#     tables_list.sort(key=_natural_keys)

#     # Order by mass
#     myorder = list(range(11, len(tables_list))) + list(range(3, 11)) + list(range(3))
#     tables_list = [tables_list[i] for i in myorder]
#     return tables_list


def load_cached_tables(orbit_sideon, cache_file='data_d.pkl', force=False):
    if os.path.isfile(cache_file) and not force:
        logger.info(f"Loaded cache {cache_file}")
        d = pickle.load(open(cache_file, 'rb'))
    else:
        tables_list = SIM_NAME_DICT.values()
        d = load_tables(tables_list, orbit_sideon)
    return d


def _get_data_dir(subdir='data'):
    return os.path.join(os.path.dirname(__file__), subdir)


def get_pickle_data(cache_file):
    fullpath = os.path.join(_get_data_dir(), cache_file)
    logger.info(f"Loading cache {fullpath}")
    d = pickle.load(open(fullpath, 'rb'))
    return d


def get_big_df(cache_file):
    return pd.concat([v for v in get_pickle_data(cache_file).values()], axis=0)


def load_tables(sim_name_list, orbit_sideon):
    """Load dataframes in a dictionary"""
    d = dict()
    for sim_name in sim_name_list:
        if orbit_sideon:
            sim_name += '_orbit_sideon'
#         print(sim_name)
        d[shorten_name(sim_name)] = get_df(sim_name)
    return d


def get_tables(sim_name, orbit_sideon, data_dir=DATA_DIR):
    appendix = "" if not orbit_sideon else "_orbit_sideon"
    filename = sim_name + appendix + '.fits'
    name = os.path.join(data_dir, "tables", filename)
    logger.debug("Getting tables: {}".format(name))
    return Table.read(name)


def get_phot(sim_name, orbit_sideon, data_dir=DATA_DIR):
    appendix = "" if not orbit_sideon else "_orbit_sideon"
    filename = sim_name + appendix + '_photometry.fits'
    logger.debug("Getting photometry: {}".format(filename))
    tbl = Table.read(os.path.join(data_dir, 'photometry', filename))
    return tbl


def get_maps(sim_name, orbit_sideon, data_dir=DATA_DIR):
    if orbit_sideon:
        tbl = Table.read(os.path.join(data_dir, 'maps_orbit_sideon_sig_los', sim_name+'_orbit_sideon_maps.fits'))
    else:
        tbl = Table.read(os.path.join(data_dir, 'maps_orbit_faceon_sig_los', sim_name+'_maps.fits'))
    return tbl


def get_maps_all_band(sim_name, orbit_sideon, data_dir=DATA_DIR):
    if orbit_sideon:
        tbl = Table.read(os.path.join(data_dir, 'maps_orbit_sideon_sig_los', sim_name+'_orbit_sideon_maps_allbands.fits'))
    else:
        tbl = Table.read(os.path.join(data_dir, 'maps_orbit_faceon_sig_los', sim_name+'_maps_allbands.fits'))
    return tbl


def get_maps_HI(sim_name, orbit_sideon, data_dir=DATA_DIR):
    appendix = "" if not orbit_sideon else "_orbit_sideon"
    filename = sim_name + appendix + '_HI_maps.fits'
    tbl = Table.read(os.path.join(data_dir, 'hi_maps', filename))
    return tbl


def get_HI_data(sim_name, orbit_sideon, data_dir=DATA_DIR):
    appendix = "" if not orbit_sideon else "_orbit_sideon"
    filename = sim_name + appendix + '_HI.fits'
    tbl = Table.read(os.path.join(data_dir, 'hi', filename))
    return tbl


def get_angmom(sim_name, orbit_sideon, data_dir=DATA_DIR):
    appendix = "" if not orbit_sideon else "_orbit_sideon"
    filename = sim_name + appendix + '_derot_angmom.fits'
    # FIXME, maybe folder not necessary
    folder = "angmom_orbit_faceon" if not orbit_sideon else 'angmom_orbit_sideon'
    tbl = Table.read(os.path.join(data_dir, folder, filename))
    return tbl


def get_magnitudes(sim_name, data_dir=DATA_DIR):
    tbl = Table.read(os.path.join(data_dir, 'magnitudes', sim_name+'_mag.fits'))
    return tbl


def get_dm(sim_name, data_dir=DATA_DIR):
    tbl = Table.read(os.path.join(data_dir, 'dm', sim_name+'_dm.fits'))
    return tbl

def get_sf(sim_name, data_dir=DATA_DIR):
    tbl = Table.read(os.path.join(data_dir, 'sf', sim_name+'_sf.fits'))
    return tbl

def get_cii(sim_name, data_dir=DATA_DIR):
    tbl = Table.read(os.path.join(data_dir, 'cii', sim_name+'_cii.fits'))
    return tbl

def get_tidal_radius(sim_name, data_dir=DATA_DIR):
    tbl = Table.read(os.path.join(data_dir, 'tidal_radius', sim_name+'_rt.fits'))
    return tbl


def get_virial(sim_name, data_dir=DATA_DIR):
    tbl = Table.read(os.path.join(data_dir, 'virial', sim_name+'_vr.fits'))
    return tbl


def get_lambda_r(sim_name, orbit_sideon, data_dir=DATA_DIR):
    appendix = "" if not orbit_sideon else "_orbit_sideon"
    filename = sim_name + appendix + '_lambda_r.fits'
    tbl = Table.read(os.path.join(data_dir, 'lambda_r', filename))
    return tbl


def get_structure(sim_name, orbit_sideon, data_dir=DATA_DIR):
    appendix = "" if not orbit_sideon else "_orbit_sideon"
    filename = sim_name + appendix + '_structure.fits'
    tbl = Table.read(os.path.join(data_dir, 'structure', filename))
    return tbl


def get_sigma(sim_name, orbit_sideon, data_dir=DATA_DIR):
    appendix = "" if not orbit_sideon else "_orbit_sideon"
    filename = sim_name + appendix + '_sigma.fits'
    tbl = Table.read(os.path.join(data_dir, 'sigma', filename))
    return tbl


def get_cold_gas(sim_name, data_dir=DATA_DIR):
    filename = sim_name + '_cold_gas.fits'
    tbl = Table.read(os.path.join(data_dir, 'cold_gas', filename))
    return tbl

# PROF_TYPES = {"dens":"dens",
#               "3d_dens": "dens", "dens"}

def get_profiles(sim_name, orbit_sideon, prof_type, data_dir=DATA_DIR):
    appendix = "" if not orbit_sideon else "_orbit_sideon"
    filename = sim_name + appendix + f'_{prof_type}.fits'
    tbl = Table.read(os.path.join(data_dir, 'profiles', filename))
    return tbl



# def _generate_sim_dict():
#     sim_d = dict()
#     for simpath in get_sim_name_list():
#         print(simpath)
#         sim_name = get_sim_name(os.path.join(simpath, 'out'))
#         sim_d[shorten_name(simpath)] = os.path.splitext(sim_name)[0]
#     return sim_d


def short2simname(name, sim_names_dict=SIM_NAME_DICT):
    return sim_names_dict.get(name)


if __name__ == '__main__':
    import pprint
    NO_GAS_DICT = dict()
    no_gas_list = glob.glob('/home/michele/sim/MySimulations/ok_new_adhoc_or_not_affected/*no_gas')
    for simpath in no_gas_list:
        print(simpath)
        sim_name = get_sim_name(os.path.join(simpath, 'out'))
        NO_GAS_DICT[shorten_name(simpath)] = sim_name
    pprint.pprint(NO_GAS_DICT)