#!/usr/bin/env python3 
import subprocess
import os
import snap_io
import argparse
import pynbody
import tqdm

parser = argparse.ArgumentParser()
parser.add_argument('snap_dir', help='Input file')
parser.add_argument('-o', '--output_dir', '--out', help='Output directory (create it if it does not exist)')
parser.add_argument('-f', '--family', help='Family of particle to use to recenter', default=None)
args = parser.parse_args()

snap_list = snap_io.load_sim(args.snap_dir)

os.makedirs(args.output_dir, exist_ok=True)

for snap in tqdm.tqdm(snap_list):
    if args.family is not None:
        snap = snap.__getattr__(args.family)
    pynbody.analysis.halo.center(snap)

    snap.write(filename=os.path.join(args.output_dir, os.path.basename(snap.filename)))
