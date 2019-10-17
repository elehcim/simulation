import os
import simulation
import pynbody
import pandas as pd
import matplotlib.pylab as plt
import matplotlib
import numpy as np
import tqdm
from astropy.table import Table
from simulation.util import make_lowess
import glob
from simulation.derived import feh, mgfe, gas_metals, neutral_fraction
import pickle
from simulation.vlos_profiles import get_max_vlos


DATA_DIR = '/home/michele/sim/analysis/ng_ana/data'
TABLE_LIST_DIR = '/home/michele/sim/analysis/ng_ana/data/tables/general/'

def get_name_peri(sim_name):
    mb, peri = sim_name.split('_')[:2]
    return mb[3:], int(peri[1:])


def get_radial_period(sim_name, data_dir=DATA_DIR):
    """Return radial period in internal code units"""
    mb, peri = sim_name.split('_')[:2]
    filename = os.path.join(data_dir, "radial_period", "radial_period_{}.pickle".format(mb))
    # print(filename)
    d = pickle.load(open(filename, 'rb'))
    # print(d)
    rperiod = d[int(peri[1:])]
    return rperiod

def compute_t_period(sim_name):
    df = get_tables(sim_name).to_pandas()
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
    df['orbital_phase'] = pd.cut(df['t_period'], [-0.25, 0.25, 0.75, 1.25, 1.75], labels=False)
    return df

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
    df = Table.read(os.path.join(data_dir, 'tables/{}.fits'.format(sim_name))).to_pandas()

    # Remove wrong sigmas
    unwanted_columns = ['sigma_star', 'sigma_gas']
    for u in unwanted_columns:
        if u in df.columns: del df[u]

    is_sideon = is_orbit_sideon(sim_name)
    name_no_orientation = sim_name.split('_orbit_sideon')[0]

    phot_tbl = Table.read(os.path.join(data_dir, 'photometry/{}_photometry.fits'.format(sim_name))).to_pandas()

    struct_tbl = Table.read(os.path.join(data_dir, "structure/{}_structure.fits".format(sim_name))).to_pandas()
    dm_tbl = Table.read(os.path.join(data_dir, "dm/{}_dm.fits".format(name_no_orientation))).to_pandas()
    mag_tbl = Table.read(os.path.join(data_dir, "magnitudes/{}_mag.fits".format(name_no_orientation))).to_pandas()
    sig_tbl = Table.read(os.path.join(data_dir, "sigma/{}_sigma.fits".format(sim_name))).to_pandas()
    cg_tbl = Table.read(os.path.join(data_dir, "../gas/cold_gas_data/{}_cold_gas.fits".format(name_no_orientation))).to_pandas()
    sf_tbl = Table.read(os.path.join(data_dir, "sf/{}_sf.fits".format(name_no_orientation))).to_pandas()

    # Merge data
    for col in phot_tbl.columns:
        df[col] = phot_tbl[col]

    for col in dm_tbl.columns:
        df[col] = dm_tbl[col]

    for col in sig_tbl.columns:
        df[col] = sig_tbl[col]

    for col in cg_tbl.columns:
        df[col] = cg_tbl[col]

    for col in mag_tbl.columns:
        df['mag_' + col] = mag_tbl[col]

    for col in ('n', 'mu_e'):
        df[col] = struct_tbl[col]

    for col in ('medians', 'averages', 'length', 'bx', 'by', 'bz'):
        df['sf_' + col] = sf_tbl[col]

    df['r_eff_fit'] = struct_tbl['r_eff']

    df['name'], df['pericenter'] = get_name_peri(sim_name)

    avg_columns = ['lambda_r', 'n', 'sfr', 'r_eff', 'r_eff3d',
                   'r_eff_fit', 'ellipticity', 'mu_e',
                   'sigma_gas', 'sigma_star',
                   'mass_star', 'orientation', 'dm_mass', 'dm_mass_all']

    # Max v_los
    try:
        df['max_v_los'] = get_max_vlos(sim_name, is_sideon)
        avg_columns.append('max_v_los')
    # FIXME ok but this is not what I need
    except Exception as e:
        print(e)
        print("Not computing Max vlos")

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
        radial_period = get_radial_period(sim_name, data_dir)
        zero_crossings = np.where(np.diff(np.signbit(df['r'].diff())))[0]
        idx_peri = zero_crossings[1]  # first element is always 0, so the second is the actual pericenter
        first_pericenter_time = df.t[idx_peri]
        # print(zero_crossings)
        # print(first_pericenter_time)

        # Put 0 of the scale on first pericenter
        df['t_period'] = (df.t-first_pericenter_time)/radial_period
        df['orbital_phase'] = pd.cut(df.t_period, [-0.25, 0.25, 0.75, 1.25, 1.75], labels=False)
            # labels=(1, 2, 3, 4))
            # labels=('peri_1', 'apo_1', 'peri_2', 'apo_2' ))
        # 0 on start of simulation
        # df['t_period'] = (df.t-df.t.iloc[0])/radial_period

    except Exception as e:
        print(e)
        print("Not computing Radial period")

    print('{}: ok, data loaded'.format(sim_name))
    return df


def get_sf_pos(sim_name):
    return pickle.load(open(os.path.join(DATA_DIR, 'sf', sim_name + "_sf_pos.pkl"), 'rb'))

def get_mach(sim_name, data_dir=DATA_DIR):
    name = os.path.join(data_dir, f"mach/{sim_name}_mach.fits")
    print(f"Getting Mach number table: {name}")
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
    return col

def get_styles_from_name(name, scatter=False):
    if scatter:
        styles = ['o', '^', "s", 'd']
    else:
        styles = ['-', '--', '-.', ':']
    if name.endswith('50'):
        st = styles[0]
    elif name.endswith('100'):
        st = styles[1]
    elif name.endswith('200'):
        st = styles[2]
    elif name.endswith('300'):
        st = styles[3]
    return st

def get_color_styles(d, scatter=False):
    col_d = dict()
    for k in d.keys():
        col_d[k] = get_color_from_name(k)
    st_d = dict()
    for k in d.keys():
        st_d[k] = get_styles_from_name(k, scatter=scatter)
    return col_d, st_d

# From https://www.tutorialspoint.com/How-to-correctly-sort-a-string-with-a-number-inside-in-Python
import re

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split('(\d+)',text) ]

def is_orbit_sideon(sim_name):
    if "_orbit_sideon" in sim_name:
        return True
    else:
        return False

def get_sim_name_list(table_list_dir=TABLE_LIST_DIR):
    tables_list = glob.glob(os.path.join(table_list_dir, '*.fits'))
    tables_list.sort(key=natural_keys)

    # Order by mass
    myorder = list(range(11, len(tables_list))) + list(range(3, 11)) + list(range(3))
    tables_list = [tables_list[i] for i in myorder]
    return tables_list

def load_cached_tables(orbit_sideon, cache_file='data_d.pkl', force=False):
    if os.path.isfile(cache_file) and not force:
        print(f"Loaded cache {cache_file}")
        d = pickle.load(open(cache_file, 'rb'))
    else:
        tables_list = glob.glob(os.path.join(TABLE_LIST_DIR, '*.fits'))
        tables_list.sort(key=natural_keys)

        # Order by mass
        myorder = list(range(11, len(tables_list))) + list(range(3, 11)) + list(range(3))
        tables_list = [tables_list[i] for i in myorder]
        d = load_tables(tables_list, orbit_sideon)
    return d

def load_tables(file_list, orbit_sideon):
    """Load dataframes in a dictionary"""
    d = dict()
    for f in file_list:
        sim_name = os.path.splitext(os.path.basename(f))[0]
        if orbit_sideon:
            sim_name += '_orbit_sideon'
#         print(sim_name)
        d[shorten_name(f)] = get_df(sim_name)
    return d


def get_tables(sim_name, data_dir=DATA_DIR):
    name = os.path.join(data_dir, "tables/{}.fits".format(sim_name))
    print("Getting tables: {}".format(name))
    return Table.read(name)


def get_phot(sim_name, data_dir=DATA_DIR):
    name = os.path.join(data_dir, "photometry/{}_photometry.fits".format(sim_name))
    print("Getting photometry: {}".format(name))
    return Table.read(name)


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
    return vlos_map