# plot COG

import argparse
import os
import collections
import numpy as np
import pynbody
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from util import first_last_snap, snapshot_list

parser = argparse.ArgumentParser()

parser.add_argument("--dir", default='~/sim/nfw_negative/out')
parser.add_argument('-n','--snap', default=None)
parser.add_argument('--npy-file', help="Specify a file with cached COG data")

args = parser.parse_args()
args.dir = os.path.expanduser(args.dir)


snaplist = snapshot_list(args.dir, include_dir=True)
print(snaplist)
first_snap, last_snap = first_last_snap(args.dir)
print("Found snapshots [{}: {}]".format(first_snap, last_snap))

# snapshots = range(first_snap, last_snap)

def compute_cog(snapshots, directory, save_cache=True):
	if not isinstance(snapshots, collections.Iterable):
		snapshots = [snapshots]

	cog = np.zeros((3, len(snapshots)), dtype=float)
	fdir = os.path.expanduser(args.dir)
	i = 0
	for snap in snapshots:
		sim = pynbody.load(os.path.join(directory, snap))
		print("reading file {}".format(snap))
		print("snapshot time {} Gyr".format(sim.properties['time'].in_units('Gyr')))

		mass = sim['mass']
		pos = sim['pos']

		tot_mass = mass.sum()
		# print "tot mass = {} Msol".format(tot_mass)
		cog[:,i] = np.array([(pos[:,0] * mass).sum(), (pos[:,1] * mass).sum(), (pos[:,2] * mass).sum()])/tot_mass

		i += 1
	print(cog)
	if save_cache:
		np.save("cog_pynbody.npy", cog)
	return cog

if __name__ == '__main__':
	if not args.npy_file:
		cog = compute_cog(snaplist, args.dir)
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

