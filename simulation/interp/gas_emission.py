import os
import sys
import argparse
import numpy as np
import pynbody
from scipy.interpolate import interpn, RegularGridInterpolator

_TABLE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'tables')

def load_data(table_type):
    with np.load(os.path.join(_TABLE_DIR, f'{table_type}_data_all.npz')) as tbl_data:
        data        = tbl_data['data']
        dens_conv   = tbl_data['dens_conv']
        temp        = tbl_data['temp']
        feh         = tbl_data['feh']
        mgfe        = tbl_data['mgfe']
        z           = tbl_data['z']
        dens        = tbl_data['dens']

    return data, dens_conv, temp, feh, mgfe, z, dens


class TableInterpolator:
    def __init__(self, table_type, bounds_error=False, fill_value=None):
        self.table_type = table_type
        self._table_data, self._dens_conv, self._temp, self._feh, self._mgfe, self._z, self._n_H = load_data(table_type)
        self._dens_conv_interpolator = RegularGridInterpolator(points=(self._feh, self._mgfe), values=self._dens_conv, bounds_error=bounds_error, fill_value=fill_value)
        self._interpolator = RegularGridInterpolator(points=(self._temp, self._feh, self._mgfe, self._z, self._n_H), values=self._table_data, bounds_error=bounds_error, fill_value=fill_value)

    # def interpolate(self, temp, feh, mgfe, z, rho):
    #     # in zero metallicity case (FeH=-99), set MgFe to some value within the available range
    #     # (specific value does not matter, since in FeH=-99 case the MgFe
    #     # value does not matter, table values are the same)

    #     # rho should be in g/cm3. My tables are in amu/cm3 (H number density), so first I have to
    #     # convert g/cm3 to amu/cm3. This conversion depends in turn on MgFe and FeH, so I first need
    #     # to interpolate on the conversion factor table toget the correct factor.

    #     if (feh == -98):
    #         mgfe = self._mgfe[0];

    #     conversion_factor = interpn((self._feh, self._mgfe), self._dens_conv, (feh, mgfe), bounds_error=False, fill_value=None)
    #     print('Density conversion factor: {}'.format(conversion_factor[0]))
    #     n_H = rho * conversion_factor[0]

    #     query = [temp, feh, mgfe, z, n_H]
    #     print('Final value to use: temp={} K, feh={}, mgfe={}, z={}, n_H (1/cm3)={}'.format(*query))

    #     result = interpn((self._temp, self._feh, self._mgfe, self._z, self._n_H), self._table_data, query, bounds_error=False, fill_value=None)
    #     return result

    def interpolate_vec(self, temp, feh, mgfe, z, rho, hyplot_boundary=False):
        # in zero metallicity case (FeH=-99), set MgFe to some value within the available range
        # (specific value does not matter, since in FeH=-99 case the MgFe
        # value does not matter, table values are the same)

        # `rho` should be in g/cm3. My tables are in amu/cm3 (H number density), so first I have to
        # convert g/cm3 to amu/cm3. This conversion depends in turn on MgFe and FeH, so I first need
        # to interpolate on the conversion factor table toget the correct factor.

        assert len(temp) == len(feh) == len(mgfe) == len(z) == len(rho)
        if isinstance(rho, pynbody.array.SimArray):
            rho_g_cm3 = rho.in_units('g cm**-3').view(np.ndarray)
        else:
            rho_g_cm3 = rho
        mymgfe = mgfe.copy()

        mymgfe[np.where(feh == -98)] = self._mgfe[0]
        # mymgfe[np.where(mgfe > self._mgfe[-1])] = self._mgfe[-1]

        mytemp = temp.copy()
        myz = z.copy()
        if hyplot_boundary:
            mytemp[np.where(temp > self._temp[-1])] = self._temp[-1]
            mytemp[np.where(temp < self._temp[0])] = self._temp[0]
            myz[np.where(z > self._z[-1])] = self._z[-1]
            myz[np.where(z < self._z[0])] = self._z[0]

        query = np.array([mytemp.view(np.ndarray),
                feh.view(np.ndarray),
                mymgfe.view(np.ndarray),
                myz.view(np.ndarray),
                rho_g_cm3])

        # print('temp={}, feh={}, mgfe={}, z={}, rho={}'.format(*query))
        # print(mymgfe.max())
        conversion_factor = self._dens_conv_interpolator(xi=(feh, mymgfe))
        # print('Density conversion factor: {}'.format(conversion_factor))
        # print('Density conversion factor shape: {}'.format(conversion_factor.shape))
        # print(query)
        query[-1] *= conversion_factor[0]

        # result = interpn((self._temp, self._feh, self._mgfe, self._z, self._n_H), self._table_data, query.T, bounds_error=True, fill_value=None)
        result = self._interpolator(xi=(query.T))

        # Extrapolation can bring negative values
        result[np.where(result < 0)] = 0

        return result

    def print_limits(self):
        print('Table limits:')
        formatting = '{:13s} - {:10.5g} {:10.5g}'
        for name, qty in zip((self.table_type, 'dens_conv', 'temp', 'feh', 'mgfe', 'z', 'n_H'), \
                         [self._table_data, self._dens_conv, self._temp, self._feh, self._mgfe, self._z, self._n_H]):
            print(formatting.format(name, qty.min(), qty.max()))
        print(formatting.format('rho (g/cm3)', self._n_H.min()/self._dens_conv.min(),
                                                  self._n_H.max()/self._dens_conv.max()))


_table_interpolator_CII = TableInterpolator('CII')
_table_interpolator_halpha = TableInterpolator('Halpha')