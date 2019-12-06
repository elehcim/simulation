import argparse
from simulation.simdata import SIM_NAME_DICT, save_tables

def parse_args(cli=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--on-orbit-plane', action='store_true', help='Put snapshot on orbit plane')
    args = parser.parse_args(cli)
    return args


def main(cli=None):
    args = parse_args(cli)
	save_tables(SIM_NAME_DICT.values(), orbit_sideon=args.om_orbit_plane)

if __name__ == '__main__':
	main()