import os
import tqdm
import gc
from pprint import pprint

import numpy as np
import pynbody
from astropy.table import Table, Column
from astropy import units as u

from simulation.util import setup_logger, get_sim_name
from simulation.save_maps import get_outname, single_snap_maps, parse_args, get_derotation_parameters, write_out_table

BANDS_AVAILABLE = ['u', 'b', 'v', 'r', 'i', 'j', 'h', 'k']

logger = setup_logger('save_maps_all_bands', logger_level='INFO')

COLUMNS_UNITS = dict()


def simulation_maps(sim_path, width, resolution,
                    side=True, face=None,
                    quat_dir=None, pivot=None, derotate=True, on_orbit_plane=False,
                    save_single_image=False, overwrite=False):

    if not os.path.isdir(sim_path):
        raise RuntimeError('Simulation path should be a directory')

    from simulation.util import snapshot_list
    snap_list = snapshot_list(sim_path, include_dir=True)

    sim_name = get_sim_name(sim_path)

    if derotate:
        quat_arr, omega_mb_arr, pivot = get_derotation_parameters(sim_name, quat_dir, pivot, snap_list)
    else:
        quat_arr = omega_mb_arr = pivot = None

    nan_arr = np.empty((resolution, resolution), dtype=np.float32)
    nan_arr[:] = np.nan

    maps_dict = dict()
    for band in BANDS_AVAILABLE:
        maps_dict['mag_{}'.format(band)] = list()
        maps_dict['lum_{}'.format(band)] = list()
        COLUMNS_UNITS['mag_{}'.format(band)] = u.mag * u.arcsec**-2
        COLUMNS_UNITS['lum_{}'.format(band)] = u.solLum * u.pc**-2

    data_out_name = get_outname('data', out_dir=sim_name, band='all', width=width, resolution=resolution)

    appendix = "" if not on_orbit_plane else "_orbit_sideon"
    out_name = sim_name + appendix + '_maps_allbands.fits'

    if os.path.isfile(out_name) and not overwrite:
        raise RuntimeError(f'File {out_name} already esist')

    for i, snap_name in enumerate(tqdm.tqdm(snap_list)):

        snap = pynbody.load(snap_name)
        time = snap.header.time
        quat = quat_arr[i, :] if quat_arr is not None else None
        omega_mb = omega_mb_arr[i, :] if omega_mb_arr is not None else None

        try:
            im = single_snap_maps(snap_name=snap_name,
                                  width=width,
                                  resolution=resolution,
                                  side=side,
                                  face=face,
                                  quat=quat,
                                  omega_mb=omega_mb,
                                  pivot=pivot,
                                  on_orbit_plane=on_orbit_plane,
                                  center_velocity=False,  # For imaging I do not need to center velocity.
                                  )

            for band in BANDS_AVAILABLE:
                lum, mag = im.sb(band)
                maps_dict['mag_{}'.format(band)].append(mag)
                maps_dict['lum_{}'.format(band)].append(lum)

            if save_single_image:
                mdict_loc = {'mag_{}'.format(band) : im.sb_mag(band) for band in BANDS_AVAILABLE}
                mdict_loc_lum = {'lum_{}'.format(band) : im.sb_lum(band) for band in BANDS_AVAILABLE}
                mdict_loc.update(mdict_loc_lum)

                mtbl_loc = Table(mdict_loc, meta={'time':time, 'time_u': str(u.kpc*u.s*u.km**-1)})
                for col_name, col in mtbl_loc.columns.items():
                    col.unit = COLUMNS_UNITS[col_name]

                mtbl_loc.write(os.path.join(data_out_name, 'maps{}_img.fits.gz'.format(snap_name[-4:])), overwrite=True)

        except ValueError as e:
            # Usually is 'Insufficient particles around center to get velocity'
            logger.error(snap_name)
            logger.error(e)

            for k, v in maps_dict.items():
                v.append(nan_arr)

        del snap
        snap_list[i]=None
        if i % 5 == 0:
            gc.collect()

    meta = dict(BANDS=''.join(BANDS_AVAILABLE), resol=resolution, width=width)
    write_out_table(maps_dict, out_name, meta=meta, column_units=COLUMNS_UNITS, overwrite=overwrite)



def main(cli=None):
    args = parse_args(cli)

    pprint(args.__dict__)
    im = None

    if args.snap_name is not None:
        im = single_snap_maps(snap_name=args.snap_name,
                              width=args.width,
                              resolution=args.resolution,
                              side=args.side,
                              face=args.face,
                              quat=np.array(args.quat.split(), dtype=np.float64),
                              pivot=np.array(args.pivot.split(), dtype=np.float64),
                              derotate=args.no_derot,
                              )
    elif args.sim_path is not None:
        simulation_maps(sim_path=args.sim_path,
                        width=args.width,
                        resolution=args.resolution,
                        side=args.side,
                        face=args.face,
                        quat_dir=args.quat_dir,
                        pivot=args.pivot,
                        derotate=args.no_derot,
                        on_orbit_plane=args.on_orbit_plane,
                        save_single_image=args.save_single_image,
                        overwrite=args.overwrite,
                        )

    return im


if __name__ == '__main__':
    main()
