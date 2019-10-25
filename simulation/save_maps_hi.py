import argparse
import functools
import os
import tqdm
import gc
from pprint import pprint
from collections import defaultdict

import numpy as np
import pynbody
from astropy.table import Table, Column
from astropy import units as u

from simulation.util import setup_logger, get_sim_name
from simulation.save_maps import single_snap_maps, parse_args, get_derotation_parameters, write_out_table
from simulation.simdata import get_center

logger = setup_logger('save_maps_hi', logger_level='INFO')

COLUMNS_UNITS = dict()


def simulation_hi_maps(sim_path, width, resolution,
                       side=True, face=None,
                       quat_dir=None, pivot=None, derotate=True, on_orbit_plane=False,
                       overwrite=False):

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

    maps_dict = defaultdict(list)
    COLUMNS_UNITS['sigma_hi'] = u.solMass * u.pc**-2

    appendix = "" if not on_orbit_plane else "_orbit_sideon"
    out_name = sim_name + appendix + '_HI_maps.fits'

    if os.path.isfile(out_name) and not overwrite:
        raise RuntimeError(f'File {out_name} already esist')

    center_list = get_center(sim_name)

    for i, snap_name in enumerate(tqdm.tqdm(snap_list)):

        snap = pynbody.load(snap_name)
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
                                  sphere_edge=None,
                                  center_velocity=False,  # For imaging I do not need to center velocity.
                                  center_pos=None,
                                  # center_pos=center_list[i],
                                  )

            maps_dict['sigma_hi'].append(im.sigma_hi())

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

    meta = dict(resol=resolution, width=width)
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
        simulation_hi_maps(sim_path=args.sim_path,
                        width=args.width,
                        resolution=args.resolution,
                        side=args.side,
                        face=args.face,
                        quat_dir=args.quat_dir,
                        pivot=args.pivot,
                        derotate=args.no_derot,
                        on_orbit_plane=args.on_orbit_plane,
                        overwrite=args.overwrite,
                        )

    return im


if __name__ == '__main__':
    main()
