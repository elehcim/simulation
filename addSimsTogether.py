import math
import os, sys, time
import getopt
import numpy as np

# different paths depending on if you run it on the cluster or on your computer
sys.path.append("/home/rpverbek/programs/hyplot/hyplot_build/src/sip")


gadgetDir = "/home/rpverbek/programs/gadget2/"
outputDir = "/home/rpverbek/programs/gadget2/output/"

import chyplot
import enums

# **** Parameters to be set by hand ***
# simulations that have to be put together
sims = [14001, 21001, 24001]
# snapshots of galaxies that have to be merged
snaps = [36, 36, 36]
# name of the new simulation
newsim = 10001
# x, y, z => distance where gas cloud should be placed (in kpc)
positions = [(0,0,0), (10, 10, 10), (10, -10, 10)]
# vx, vy, vz => velocities of the gas cloud (in km/sec)
velocities = [(0,0,0), (0, 2.5, 0), (-5, 3, 0)]
# *************************************

kpc = 3.08568025e19 # conversion of kpc to m
vel = 1000 # conversion of km/s to m/s
gyr = 31556926e9 # conversion of gigayears to sec
solMass = 1.98892e36 # conversion of 10**6 solar mass to kg

totalData = None
reader = chyplot.CDataGadget()
writer = chyplot.CWriteGadget()


for sim, snap, position, velocity in zip(sims, snaps, positions, velocities):
    reader.setPrefix(os.path.join(outputDir, "sim%04.d"%sim))
    reader.setRunNumber(sim)
    reader.set_file(snap)

    try:
        data = reader.readFile()
    except chyplot.IOError as e:
        print
        print "****Error reading file****"
        print "simulation number:", sim
        print e
        print e.what()
        sys.exit(12)

    data.rcom(True, enums.T_star)
    data.vcom(True, enums.T_star)

    # set the distance of the gas cloud from the galaxy: x,y,z
    x, y, z = position
    data.translate(enums.T_all,x, y, z) 
    # set the speed of the gas cloud: vx, vy, vz
    vx, vy, vz = velocity
    data.kick(enums.T_all, vx, vy, vz) 

    print "used sim", sim, " snapshot ", snap, " which has time ", data.time()

    if totalData == None:
        totalData = data
    else:
        # add datablocks together
        totalData = totalData + data


try:
    writer.writeFile(totalData, os.path.join(gadgetDir, "ICs", "sim%04.d"%newsim + ".gic"), enums.T_all)
    print "Generated new IC file, sim", newsim, "with time", totalData.time()
except chyplot.IOError as e:
    print
    print "****Error writing file****"
    print "simulation number:", newsim
    print e
    print e.what()
    sys.exit(13)


## This code generates a .gic file from two other .gic files now to do:
## change parameter file
## make directories on node
## run as: bash gogogadget.bash 5029 8 2
## don't forget the final 2 this is a restart flagg, otherwise he put all birthyears of the stars at the start of the simulation (7 Gyr for example)

