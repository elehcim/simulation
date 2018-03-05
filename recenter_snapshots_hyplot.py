# hyplot visual script to read in snapshot, recenter, and write recentered snapshot again
import sys
import subprocess
import os

import chyplot
import enums

import argparse

from util import first_last_snap


# args = sys.argv[1:] # getting relevant arguments (first two are 'hyplot', '--visual=thisScript.py')
# fdir = "/home/michele/sim/MySimulations/Moria8Gyr_tidal/sim69002_p200.0_a600.0_r600.0_c8.15_movie/out"


parser = argparse.ArgumentParser()
parser.add_argument('snap_dir', help='Snapshots directory')
parser.add_argument('-o', '--output_dir', '--out', help='Output directory (create it if it does not exist)')
parser.add_argument('-s', '--start', default=None, type=int, help='First snaphot to consider')
parser.add_argument('-e', '--end',   default=None, type=int, help='Last snaphot to consider ')
args = parser.parse_args()

# run = int(args[0])
# snapStart = int(args[1])
# snapEnd = int(args[2])
snapName = 'snapshot_'

snapStart, snapEnd = first_last_snap(args.snap_dir)

if args.start is not None:
    snapStart = args.start
if args.end is not None:
    snapEnd = args.end

print snapStart, snapEnd, snapName

dr = chyplot.CDataGadget()
dr.setPrefix( args.snap_dir )
dr.checkFilesPresent() # set the first and last dump

for snap in range(snapStart, snapEnd+1):
    dr.set_file(dr.firstDump() + snap)
    data = dr.readFile()

    print " * recentering on stars "
    data.rcom(True, enums.T_star,0, 0, 0, True)
    data.vcom(True, enums.T_star)

    chyplot.cglobals.plmap.setDataBlock(data)

    newDir = output_dir
    # newDir = "/home/michele/sim/MySimulations/Moria8Gyr_tidal/sim69002_p200.0_a600.0_r600.0_c8.15_movie_recentered"
    if not os.path.exists(newDir):
        os.system('mkdir -p %s' %newDir)
    fileName = os.path.join(newDir, snapName + "%04.d" % snap)
    print fileName
    
    writer = chyplot.CWriteGadget()
    writer.writeFile(data, fileName, enums.T_all)
