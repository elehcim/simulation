# plot COG

import argparse
import enums
import os
import sys
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import enums

parser = argparse.ArgumentParser()

# parser.add_argument('--newsim', default=0, help='A new simulation will be created. This is the name of the new simulation.')
parser.add_argument("--sim")
parser.add_argument('--dir', default='~/sim/nfw_c3/out')
parser.add_argument('-n','--snap', default=None)
args = parser.parse_args()

snapshots = range(4,99)

import chyplot
getProp = chyplot.cglobals.plmap.getSecond

cog = np.zeros((3, len(snapshots)), dtype=float)
i = 0 
for snap in snapshots:
	dr = chyplot.CDataGadget(snap)
	fdir = os.path.expanduser(args.dir)
	dr.setPrefix( fdir )

	# if args.snap is None:
	# 	dr.set_file(dr.lastDump())
	# else:
	# 	dr.set_file(args.snap)

	print "reading file ", dr.filename()
	data = dr.readFile()
	print "snapshot time {} Gyr".format(data.time())
	# data.rcom(True, enums.T_star, 0, 0, 0, True)
	# data.vcom(True, enums.T_star)
	# data.rotate(enums.T_gas, 2, True)
	# data.convertUnits()
	mass = np.array(data.getDataArray(enums.T_all, getProp('mass'), True))
	# time.append(data.header) = np.array(data.getDataArray(enums.T_all, getProp('time'), True))

	x = np.array(data.getDataArray(enums.T_all, getProp('x'), True))
	y = np.array(data.getDataArray(enums.T_all, getProp('y'), True))
	z = np.array(data.getDataArray(enums.T_all, getProp('z'), True))

	tot_mass = mass.sum()
	print "tot mass = {} 10^10 Msol".format(tot_mass)

	cog[:,i] = ((x*mass).sum()/tot_mass,
				(y*mass).sum()/tot_mass, 
				(z*mass).sum()/tot_mass)

	i += 1
	del dr  # closing file

np.save("cog.npy", cog)

fig = plt.figure()
ax = fig.add_subplot(121, projection='3d')
ax.scatter(*cog)

ax2 = fig.add_subplot(122)
ax2.scatter(*cog[:2])

plt.tight_layout()
plt.show()