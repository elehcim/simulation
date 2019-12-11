import os
import argparse
from simulation.simdata import SIM_NAME_DICT, save_tables, _get_data_dir, _pickle_filename


def parse_args(cli=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--no-on-orbit-plane', action='store_false', help='Put snapshot on orbit plane')
    parser.add_argument('--no-gas', action='store_true', help='Get the no gas (only tides) version')
    args = parser.parse_args(cli)
    return args


def main(cli=None):
    args = parse_args(cli)
    orbit_sideon = args.no_on_orbit_plane
    cache_file = _pickle_filename(orbit_sideon, args.no_gas)
    fullpath = os.path.join(_get_data_dir(), cache_file)
    print(f'Saving data to {fullpath}')
    save_tables(SIM_NAME_DICT.values(), orbit_sideon=orbit_sideon)

if __name__ == '__main__':
    main()