#!/usr/bin/env python3 
import subprocess
import os
import snap_io
import argparse
import pynbody
import tqdm

parser = argparse.ArgumentParser()
parser.add_argument('snap_dir', help='Snapshots directory')
parser.add_argument('-o', '--output_dir', '--out', help='Output directory (create it if it does not exist)')
parser.add_argument('-f', '--family', help='Family of particle to use to recenter', default=None)
# parser.add_argument('-s', '--start', default=None, type=int, help='First snaphot to consider')
# parser.add_argument('-e', '--end',   default=None, type=int, help='Last snaphot to consider ')
args = parser.parse_args()

snap_list = snap_io.load_sim(args.snap_dir)

# if args.start is not None:
#     snap_list = args.start
# if args.end is not None:
#     snapEnd = args.end

os.makedirs(args.output_dir, exist_ok=True)
for snap in tqdm.tqdm(snap_list):
    if args.family is not None:
        snap = snap.__getattr__(args.family)
    pynbody.analysis.halo.center(snap)

	#FIXME Test if output file is readable by gent code
    snap.write(filename=os.path.join(args.output_dir, os.path.basename(snap.filename)))
