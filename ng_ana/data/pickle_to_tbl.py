from dump_features import make_df, make_astropy_table
from astropy.table import Table
import pickle
import astropy.units as u
import glob


if __name__ == '__main__':

    pkl = sorted(glob.glob('../moria*.pickle'))

    # my_pkl = 'mb.69002_p200_a800_r600_s10.pickle'
    # tbl = make_astropy_table(my_pkl)

    for pickle_file in pkl:
        tbl = make_astropy_table(pickle_file)
        tbl.write(pickle_file.replace('pickle', 'fits'), format='fits')
