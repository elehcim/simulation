# plot COG

import argparse
import os
import sys
import collections
import numpy as np
import astropy.units as u
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from util import first_last_snap
import enums

parser = argparse.ArgumentParser()

parser.add_argument("--dir", default='~/sim/nfw_at_the_end/out')
parser.add_argument('-n','--snap', default=None)
parser.add_argument('--npy-file', help="Specify a file with cached COG data")

args = parser.parse_args()
args.dir = os.path.expanduser(args.dir)

import chyplot
getProp = chyplot.cglobals.plmap.getSecond

first_snap, last_snap = first_last_snap(args.dir)
print("Found snapshots [{}: {}]".format(first_snap, last_snap))


snapshots = range(first_snap, last_snap)

def compute_cog(snapshots, directory, save_cache=True):
	if not isinstance(snapshots, collections.Iterable):
		snapshots = [snapshots]

	cog = np.zeros((3, len(snapshots)), dtype=float)
	fdir = os.path.expanduser(args.dir)
	i = 0
	for snap in snapshots:
		dr = chyplot.CDataGadget(snap)
		dr.setPrefix( args.dir )

		# if args.snap is None:
		# 	dr.set_file(dr.lastDump())
		# else:
		# 	dr.set_file(args.snap)

		print("reading file ", dr.filename())
		data = dr.readFile()
		print("snapshot time {} Gyr".format(data.time()))
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
		# print "tot mass = {} 10^10 Msol".format(tot_mass)

		cog[:,i] = ((x*mass).sum()/tot_mass,
					(y*mass).sum()/tot_mass,
					(z*mass).sum()/tot_mass)

		i += 1
		del dr  # closing file
	print(cog)
	if save_cache:
		np.save("cog.npy", cog)
	return cog

if __name__ == '__main__':
	if not args.npy_file:
		cog = compute_cog(snapshots, args.dir)
	else:
		cog = np.load(args.npy_file)

	fig = plt.figure()
	ax = fig.add_subplot(121, projection='3d')
	ax.set_xlabel("x (kpc)")
	ax.set_ylabel("y (kpc)")
	ax.set_zlabel("z (kpc)")
	ax.scatter(*cog)

	ax2 = fig.add_subplot(122)
	ax2.set_xlabel("x (kpc)")
	ax2.set_ylabel("y (kpc)")
	ax2.scatter(*cog[:2])

	plt.axis('equal')

	plt.tight_layout()
	plt.show()

