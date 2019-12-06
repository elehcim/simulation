import argparse
from simulation.simdata import SIM_NAME_DICT, save_tables

def parse_args(cli=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--no-on-orbit-plane', action='store_false', help='Put snapshot on orbit plane')
    args = parser.parse_args(cli)
    return args


def main(cli=None):
    args = parse_args(cli)
    orbit_sideon = args.no_on_orbit_plane
    save_tables(SIM_NAME_DICT.values(), orbit_sideon=orbit_sideon)

if __name__ == '__main__':
    main()