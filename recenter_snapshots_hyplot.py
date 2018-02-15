# hyplot visual script to read in snapshot, recenter, and write recentered snapshot again
# FIXME remove hardcoded paths and use argparse 
import sys
import subprocess
import os

import chyplot
import enums

import argparse

args = sys.argv[1:] # getting relevant arguments (first two are 'hyplot', '--visual=thisScript.py')


run = int(args[0])
snapStart = int(args[1])
snapEnd = int(args[2])
snapName = 'snapshot_'

print run, snapStart, snapEnd, snapName

dr = chyplot.CDataGadget(run)
fdir = "/home/michele/sim/MoRIA/sim%05.d/"%run
dr.setPrefix( fdir )
dr.checkFilesPresent() # set the first and last dump

for snap in range(snapStart, snapEnd+1):
    dr.set_file(dr.firstDump() + snap)
    data = dr.readFile()

    print " * recentering on stars "
    data.rcom(True, enums.T_star,0, 0, 0, True)
    data.vcom(True, enums.T_star)

    chyplot.cglobals.plmap.setDataBlock(data)

    newDir = "/home/michele/sim/recentered/sim" + str(run) +"/snaps/"
    if not os.path.exists(newDir):
        os.system('mkdir -p %s' %newDir)
    fileName = newDir + snapName + "%04.d" %snap
    print fileName
    
    writer = chyplot.CWriteGadget()
    writer.writeFile(data, fileName, enums.T_all)
