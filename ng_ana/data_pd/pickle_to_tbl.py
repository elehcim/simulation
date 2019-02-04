from dump_features import make_df
from astropy.table import Table
import pickle
import astropy.units as u
import glob
def make_astropy_table(pickle_file):
    times, mass_star, sigma_star, sigma_gas, r_eff, sfr = pickle.load(open(pickle_file, 'rb'))
    tbl = Table({'t': times * u.Gyr,
                 'mass_star':mass_star * u.solMass,
                 'sigma_star':sigma_star * u.km/u.s,
                 'sigma_gas':sigma_gas * u.km/u.s,
                 'r_eff':r_eff * u.kpc,
                 'sfr':sfr * u.solMass/u.yr})
    return tbl


if __name__ == '__main__':

    pkl = sorted(glob.glob('*.pickle'))

    # my_pkl = 'mb.69002_p200_a800_r600_s10.pickle'
    # tbl = make_astropy_table(my_pkl)

    for pickle_file in pkl:
        tbl = make_astropy_table(pickle_file)
        tbl.write(pickle_file.replace('pickle', 'fits'), format='fits')
