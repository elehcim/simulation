#!/bin/bash

SNAP_DIR=/home/michele/sim/MySimulations/hi_osc/mb.69002_p200_a800_r600/out
OUT_DIR=ssam.69p2.xy.sph

mkdir -p $OUT_DIR


# cat in | parallel -j4 -I% python oo_lambda.py --snap=% --side --out-dir=$OUT_DIR
ls $SNAP_DIR/snapshot_0??[02468] | parallel -k -j4 python oo_lambda.py --snap={} --out-dir=$OUT_DIR 2> $OUT_DIR/log.log > $OUT_DIR/tl.dat

# for i in `ls $SNAP_DIR/snapshot_???[05]`; do
# 	python oo_lambda.py --snap=$i --side --out-dir=$OUT_DIR
# done
