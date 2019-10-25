import os
import pandas as pd
import numpy as np
from astropy.table import Table

DATA_FOLDER = os.path.join(os.path.dirname(__file__), 'observationalData/')

Y_P = 0.2551  # Primordial helium abundance (Izotov 2014)
# (M_HI + M_He) / M_HI = 1.342 = 1/(1-Y_P)

def load_papastergis16(data_folder=DATA_FOLDER):
    filepath = os.path.join(data_folder, 'papastergis.txt')
    # df = pd.read_csv(filepath,
    #             comment='#',
    #             names="name,sample,Rout,Vout,Vout_err_low,Vout_err_high,Vrot,Mstar,Mbar,Vh_nfw,Mh_nfw,Vh_dc14,Mh_dc14".split(','),
    #             index_col='name',
    # #             header=4,
    # #             skiprows=7
    #            )
    # df['mHI'] = np.log10((10**df.Mbar-10**df.Mstar)/1.3)

    names = "name,sample,Rout,Vout,Vout_err_low,Vout_err_high,Vrot,Mstar,Mbar,Vh_nfw,Mh_nfw,Vh_dc14,Mh_dc14".split(',')
    units = ('', '','kpc', 'km/s','km/s','km/s','km/s','log10(Msun)','log10(Msun)','km/s','log10(Msun)','km/s','log10(Msun)')
    names_units = dict(zip(names, units))
    tbl = Table.read(filepath,
                     format='csv',
                     comment='#',
                     names=names,
                     # index_col='name',
                    )
    for cn in tbl.colnames:
        tbl[cn].unit = names_units[cn]

    tbl['mHI'] = np.log10((10**tbl['Mbar']-10**tbl['Mstar'])/1.342)
    tbl['mHI'].unit = 'log10(Msun)'

    return tbl
