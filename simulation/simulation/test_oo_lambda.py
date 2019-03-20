from simulation.oo_lambda import main
from astropy.table import Table
import astropy.units as u

ssam = main(['--snap', '~/sim/MySimulations/ng/mb.69002_p200_a800_r600/out/snapshot_0043'])


# tbl = Table(([ssam.sb_mag * u.mag/u.arcsec**2, ssam.v_los_map * u.km/u.s], [ssam.sb_mag * u.mag/u.arcsec**2, ssam.v_los_map * u.km/u.s]))
# tbl.write('asads.fits', overwrite=True)


import numpy as np
from astropy.table import Table, Column

tabname_before = '../../_maps_69p2_2/maps_data_v_w10_r400_a20.fits_img.fits.gz'

# tbl_after = transpose_table(tbl_before)