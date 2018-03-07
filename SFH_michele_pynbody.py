import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import math
import matplotlib
import argparse
import pynbody


from util import first_last_snap, get_snapshot_data

use_cosmolopy = True if sys.version_info[0] < 3 else False

matplotlib.rc('font', size=22)

maxR = 3  # Maximum radius within which you select stars

timebinslength = 0.05 # Time resolution in Gyr with which you want the star formation history

parser = argparse.ArgumentParser()
parser.add_argument("--sim", dest="simulations", default=('~/sim/sim60003',), nargs='+')
parser.add_argument('-n','--snap', default=None, type=int)
parser.add_argument('--max-time', default=None, type=float)
parser.add_argument('--no-pop3', action='store_true')
parser.add_argument('--sfr', action='store_true')

args = parser.parse_args()

timeMax = args.max_time # 13.15 # Maximum time

fig = plt.figure(figsize=(15, 8))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

axs = [ax1, ax2]
for simulation in args.simulations:

	# We are getting the SFH of all the stars in the simulation and then only within maxR
	maxRadiuses = [100000, maxR] 

	# Read the data
	fdir = os.path.expanduser(simulation)
	if args.snap is None:
		first_snap, last_snap = first_last_snap(fdir)
		print("Found snapshots [{}: {}]".format(first_snap, last_snap))
		my_snap = last_snap
	else:
		my_snap = args.snap

	sim = pynbody.load(os.path.join(fdir, "snapshot_{:04d}".format(my_snap)))
	sim.physical_units()

	pynbody.analysis.halo.center(sim)

	if timeMax is None:
		print("snapshot time {} Gyr".format(sim.properties['time'].in_units('Gyr')))
		timeMax = sim.properties['time'].in_units('Gyr')

	time = np.arange(0, timeMax, timebinslength)


	# # Exclude Pop3 particles
	# if args.no_pop3:
	# 	metallicity_visitor = getProp("[Fe/H]")
	# 	n_star = len(data.getDataArray(enums.T_star, getProp('[Fe/H]'), True))
	# 	print("Number of star particles with POP3:    {}".format(n_star))
		
	# 	# Cut out POP3 stars
	# 	data.applyLimits(metallicity_visitor, -5, 100, enums.T_star)
		
	# 	n_star_without_pop3 = len(data.getDataArray(enums.T_star, getProp('[Fe/H]'), True))
	# 	print("Number of star particles without POP3: {}".format(n_star_without_pop3))
	
	# visitor = getProp("dist_to_z")

	print('sim: {}'.format(simulation))
	print('Maximum lookback time {:.2f} Gyr'.format(timeMax))
	for maxRadius, ax in zip(maxRadiuses, axs):

		filt = sim.s['rxy'].in_units('kpc') < maxRadius
		# print(filt)
		birthtimes = sim.s[filt]['tform']
		print("{:3d} star particles within radius {:g} kpc".format(len(birthtimes), maxRadius))
		masses = sim.s[filt]['massform']
		Mstar = sum(masses)

		SFR = np.zeros(len(time))

		# Determine the number of stars being formed in each time bin, and so the star formation rate
		for birthtime, mass in zip(birthtimes, masses):
			for idx in range(len(time)):
				if time[idx] - birthtime < timebinslength and time[idx] - birthtime > 0:
					SFR[idx] += mass*10**6 / (timebinslength*10**9)


		# get the cumulative star formation history
		if args.sfr:
			SFRcum = SFR
		else:
			SFRcum = np.cumsum(SFR)

		SFRcum = SFRcum * timebinslength * 10**3 / Mstar
		# SFRcum = [s*(timebinslength*10**3)/Mstar for s in SFRcum]


		lookbacktime = [timeMax - t for t in time]
		print('Mstar = {:.4f}  at maxRadius = {}'.format(Mstar, maxRadius))

		ax.plot(lookbacktime, SFRcum, 'b-')

		ax.set_xlabel('$t_\mathrm{lookback}\ \mathrm{[Gyr]}$')
		ax.set_ylabel('$f_\star(R<{}~\mathrm{{kpc}})$'.format(maxRadius))

		ax.set_xlim(timeMax, 0)
		ax.set_ylim(0,1.)


# Plot the case of a constant SFR over time
ax1.plot([timeMax, 0], [0, 1], '--', color='k', linewidth=2)
ax2.plot([timeMax, 0], [0, 1], '--', color='k', linewidth=2)

fs = 32
ax1.text(0., 1.13, '$\mathbf{a.}$', color='k', fontsize=fs, transform=ax1.transAxes)
ax2.text(0., 1.13, '$\mathbf{b.}$', color='k', fontsize=fs, transform=ax2.transAxes)

if use_cosmolopy:
	from cosmolopy import distance
	from cosmolopy import parameters as cp
	from cosmolopy import constants as cc
	#Add a redshift axis on top
	redshifts = [13, 6, 3, 2, 1, 0]
	cosmology = cp.WMAP7_BAO_H0_mean(flat=True)

	Age0 = distance.age_flat(13.5, **cosmology)

	times = [(distance.age_flat(z, **cosmology) - Age0)/cc.Gyr_s for z in redshifts]

	ax4 = ax1.twiny()
	ax4.set_xticks(times)
	ax4.set_xticklabels(['${}$'.format(z) if not z==redshifts[0] else '' for z in redshifts])
	ax4.set_xlabel('$z$')


	ax5 = ax2.twiny()
	ax5.set_xticks(times)
	ax5.set_xticklabels(['${}$'.format(z) if not z==redshifts[0] else '' for z in redshifts])
	ax5.set_xlabel('$z$')
else:
	from pynbody.analysis import pkdgrav_cosmo as cosmo
	# Use default cosmology
	# do not use sim as input because sim.properties['omegaM0'] is 1.0 and it should be 0.272
	c = cosmo.Cosmology()

	ax4 = ax1.twiny()
	ax5 = ax2.twiny()
	for ax in (ax4, ax5):
		labelzs = [6, 3, 2, 1, 0] #np.arange(5, int(sim.properties['z']) - 1, -1)
		times = [13.5 * c.Exp2Time(1.0 / (1 + z)) / c.Exp2Time(1) for z in labelzs]
		ax.set_xticks(times)
		ax.set_xticklabels([str(x) for x in labelzs])
		ax.set_xlabel('$z$')


fig.subplots_adjust(left = 0.1, right = 0.98, bottom = 0.12, top = 0.85, hspace=0, wspace=0.25)

directory = '/home/rpverbek/programs/results/SFR/'
name = 'figureSFH_test.pdf'


plt.show()
# plt.savefig(name=directory+name, dpi=300, show=False)
# finalize(name=directory+name, dpi=300, show=False)
