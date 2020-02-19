import os
import sys
import argparse
import numpy as np
from scipy.interpolate import interpn


def load_data(filename):
    with np.load(filename) as data:
        table_data  = data['table_data']
        dens_conv   = data['dens_conv']
        temp        = data['temp']
        feh         = data['feh']
        mgfe        = data['mgfe']
        z           = data['z']
        dens        = data['dens']

    return table_data, dens_conv, temp, feh, mgfe, z, dens


# Emission is a function of temp, FeH, MgFe, Z, dens
_table_data, _dens_conv, _temp, _feh, _mgfe, _z, _dens = load_data(os.path.dirname(__file__) +'/tables/HI_data_all.npz')


def print_limits():
    print('Table limits:')
    for name, qty in zip(('HI', 'dens_conv', 'temp', 'feh', 'mgfe', 'z', 'dens'), [_table_data, _dens_conv, _temp, _feh, _mgfe, _z, _dens]):
        print('{:10s} - {:10.5g} {:10.5g}'.format(name, qty.min(), qty.max()))


def example_interpolation():
    conversion_factor = interpn((_feh, _mgfe), _dens_conv, (0, 0), bounds_error=False, fill_value=None)
    query = [10.0, 0, 0, 0, 0.01]
    assert conversion_factor[0]==4.69599000e+23
    print('conversion_factor', conversion_factor)

    query[-1] *= conversion_factor[0]

    print(query)

    result = interpn((_temp, _feh, _mgfe, _z, _dens), _table_data, query, bounds_error=False, fill_value=None)
    print("HI: {}".format(result))
    # print(result[0]-8.16961943e-40)
    # np.testing.assert_array_almost_equal_nulp(result[0], 8.16961943e-40, nulp=10)
    # assert result[0]==8.16961943e-40


def get_HI_vec(temp, feh, mgfe, z, rho):
    # in zero metallicity case (FeH=-99), set MgFe to some value within the available range
    # (specific value does not matter, since in FeH=-99 case the MgFe
    # value does not matter, table values are the same)
    mymgfe = mgfe.copy()
    mymgfe[np.where(feh == -98.0)] = _mgfe[0]

    query = np.array([temp.in_units('K').view(np.ndarray),
            feh.view(np.ndarray),
            mymgfe.view(np.ndarray),
            z*np.ones_like(feh).view(np.ndarray),
            rho.in_units('g cm**-3').view(np.ndarray)])
    # print('temp={}, feh={}, mgfe={}, z={}, rho={}'.format(*query))

    conversion_factor = interpn((_feh, _mgfe), _dens_conv, (feh, mymgfe), bounds_error=False, fill_value=None)
    # print('Density conversion factor: {}'.format(conversion_factor))
    # print('Density conversion factor shape: {}'.format(conversion_factor.shape))
    # print(query[-1])
    query[-1] *= conversion_factor[0]
    # Adjust the density
    query[-1][np.where(query[-1] > _dens.max())] = _dens.max()

    # if query[-1] > _dens.max():
    #     query[-1] = _dens.max()
    # print(query)

    result = interpn((_temp, _feh, _mgfe, _z, _dens), _table_data, query.T, bounds_error=False, fill_value=None)
    return result


if __name__ == '__main__':
    print_limits()
    parser = argparse.ArgumentParser("Get HI emission from gas particles")
    parser.add_argument('values', type=float, nargs=5, help='Five floats: Temp, FeH, MgFe, Z, density (in g/cm**3 ??)')
    args = parser.parse_args() #"10.0 0 0 0 0.01".split())

    query = args.values
    # print("Asked for:")
    # print('temp={}, feh={}, mgfe={}, z={}, rho={}'.format(*query))

    result = get_HI_vec(*query)
    # print("HI: {}".format(result[0]))
    print(result[0])