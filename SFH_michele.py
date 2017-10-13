import os
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import math
import matplotlib
import argparse

from hyplot.plot import PFigure
from hyplot.visual.PRunData import PRunData, PFileData
from util import first_last_snap
import chyplot
import enums

from cosmolopy import distance
from cosmolopy import parameters as cp
from cosmolopy import constants as cc

getProp = chyplot.cglobals.plmap.getSecond

matplotlib.rc('font', size=22)


maxR = 3  # Maximum radius within which you select stars

timebinslength = 0.05 # Time resolution with which you want the star formation history

parser = argparse.ArgumentParser()
parser.add_argument("--sim", dest="simulations", default='~/sim', nargs='+')
# parser.add_argument('--dir', default='~/sim')
parser.add_argument('-n','--snap', default=None, type=int)
parser.add_argument('--max-time', default=None, type=float)
parser.add_argument('--no-pop3', action='store_true')

args = parser.parse_args()

timeMax = args.max_time # 13.15 # Maximum time

fig = plt.figure(FigureClass = PFigure.PFigure, figsize=(15, 8))
ax1 = fig.add_my_subplot(121)
ax2 = fig.add_my_subplot(122)

axs = [ax1, ax2]

for simulation in args.simulations:

	# We are getting the SFH of all the stars in the simulation and then only within maxR
	maxRadiuses = [100000, maxR] 


	# Read the data
	dr = chyplot.CDataGadget()
	fdir = os.path.expanduser(simulation)
	dr.setPrefix( fdir )
	# dr.checkFilesPresent() # set the first and last dump # IT DOES NOT WORK!

	if args.snap is None:
		first_snap, last_snap = first_last_snap(fdir)
		print "Found snapshots [{}: {}]".format(first_snap, last_snap)
		# dr.set_file( dr.lastDump())
		dr.set_file(last_snap)
	else:
		dr.set_file(args.snap)
	
	print "reading file ", dr.filename()
	data = dr.readFile()

	data.rcom(True, enums.T_star, 0, 0, 0, True)
	data.vcom(True, enums.T_star)
	data.rotate(enums.T_gas, 2, True)

	data.convertUnits()

	if timeMax is None:
		print "snapshot time {} Gyr".format(data.time())
		timeMax = data.time()

	time = np.arange(0, timeMax, timebinslength)

	if args.no_pop3:
		# Exclude Pop3 particles
		datacopy = data.limitsCopy(visitor, 0, maxRadius, enums.T_star)
		len("Number of star particles with POP3:    ", len(datacopy))
		visitor2 = chyplot.cglobals.plmap.getSecond("[Fe/H]")
		data.applyLimits(visitor2, -5, 100, enums.T_star)
		datacopy = data.limitsCopy(visitor, 0, maxRadius, enums.T_star)
		len("Number of star particles without POP3: ", len(datacopy))


	visitor = chyplot.cglobals.plmap.getSecond("dist_to_z")

	print 'sim: {}'.format(simulation)
	print 'Maximum lookback time {:.2f} Gyr'.format(timeMax)
	for maxRadius, ax in zip(maxRadiuses, axs):

		datacopy = data.limitsCopy(visitor, 0, maxRadius, enums.T_star)

		birthtimes = datacopy.getDataArray(enums.T_star, chyplot.cglobals.plmap.getSecond('birthtime'), True)
		print "Number of star particles within radius {:8g} kpc: {}".format(maxRadius, len(birthtimes))
		masses = datacopy.getDataArray(enums.T_star, chyplot.cglobals.plmap.getSecond('initialMass'), True)

		Mstar = sum(masses)

		SFR = np.zeros(len(time))

		# Determine the number of stars being formed in each time bin, and so the star formation rate
		for birthtime, mass in zip(birthtimes, masses):
			for idx in range(len(time)):
				if time[idx] - birthtime < timebinslength and time[idx] - birthtime > 0:
					SFR[idx] += mass*10**6 / (timebinslength*10**9)


		# get the cumulative star formation history
		SFRcum = np.cumsum(SFR)
		SFRcum = [s*(timebinslength*10**3)/Mstar for s in SFRcum]


		lookbacktime = [timeMax - t for t in time]
		print 'Mstar = {:.4f}  at maxRadius = {}'.format(Mstar, maxRadius)

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


fig.subplots_adjust(left = 0.1, right = 0.98, bottom = 0.12, top = 0.85,hspace = 0, wspace=0.25)

directory = '/home/rpverbek/programs/results/SFR/'
name = 'figureSFH_test.pdf'


plt.show()
# plt.savefig(name=directory+name, dpi=300, show=False)
# finalize(name=directory+name, dpi=300, show=False)


