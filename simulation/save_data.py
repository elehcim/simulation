import os
import argparse
from simulation.simdata import SIM_NAME_DICT, NO_GAS_DICT, _get_data_dir, load_tables
import pickle
from simulation.util import setup_logger

logger = setup_logger('save_data', logger_level='INFO')

def _pickle_filename(orbit_sideon, no_gas):
    import datetime
    today = datetime.date.today().strftime("%Y%m%d")
    stem = 'data'
    if no_gas:
        stem += '_dng'
    else:
        stem += '_d'
    if orbit_sideon:
        stem += '_orbit_sideon'
    cache_file = f'{stem}_{today}.pkl'
    return cache_file


def save_tables(sim_name_list, orbit_sideon, no_gas, cache_file=None):
    """Save dictionary of dataframes in a pickle file"""
    if cache_file is None:
        cache_file = _pickle_filename(orbit_sideon, no_gas)
    fullpath = os.path.join(_get_data_dir(), cache_file)
    d = load_tables(sim_name_list=sim_name_list, orbit_sideon=orbit_sideon)
    logger.info(f'Writing data to {fullpath}')
    pickle.dump(d, open(fullpath, 'wb'))


def parse_args(cli=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--orbit-faceon', action='store_true', help='See snapshots from above the orbit plane')
    parser.add_argument('--no-gas', action='store_true', help='Get the no gas (only tides) version')
    args = parser.parse_args(cli)
    return args


def main(cli=None):
    args = parse_args(cli)
    orbit_sideon = not args.orbit_faceon
    cache_file = _pickle_filename(orbit_sideon, args.no_gas)
    fullpath = os.path.join(_get_data_dir(), cache_file)
    logger.info(f'Saving data to {fullpath}')
    if args.no_gas:
        save_tables(NO_GAS_DICT.values(), orbit_sideon=orbit_sideon, no_gas=args.no_gas)
    else:
        save_tables(SIM_NAME_DICT.values(), orbit_sideon=orbit_sideon, no_gas=args.no_gas)

if __name__ == '__main__':
    main()