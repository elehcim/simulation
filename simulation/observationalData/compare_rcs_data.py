from astropy.io import fits
import math
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator#, LogLocator
from hyplot.visual.PRunData import PFileData
import chyplot
import enums

sims = []#[14001, 41001, 57001, 60001, 62001, 67001, 68001, 69001, 71001, 72001]


fig, ax = plt.subplots(1, 1, figsize=(10,6))

inputfile = "/home/rpverbek/observationalData/papastergis.txt"

Routs = []
Vouts = []
with open(inputfile) as papastergis:
	for line in papastergis:
		if not line.startswith('#'):
			line = line.rstrip("\n")
			temp = line.split(',')
			Routs.append(temp[2])
			Vouts.append(temp[3])

ax.plot(Routs, Vouts, 'bo')


for sim in sims:

	propsVcs = PFileData('/home/rpverbek/programs/results/rotation_curves/vcs_sim{}.dat'.format(int(sim)))
	radii = propsVcs.getData('R')
	#vrots = propsVcs.getData('v_rot')
	#vas = propsVcs.getData('v_a')
	vcircs = propsVcs.getData('v_circ')

	dr = chyplot.CDataGadget(int(sim))
	fdir = "/media/DATA/simulations/sim"+str(sim)
	dr.setPrefix( fdir )
	dr.checkFilesPresent() # set the first and last dump
	dr.set_file(dr.lastDump())
	data = dr.readFile()
	data.rcom(True, enums.T_star, 0, 0, 0, True)
	data.vcom(True)
	data.convertUnits()


	galaxy = chyplot.CGalaxy(data)
	nbins = 2000
	Rmax = 50.
	vc_theors = galaxy.getCircularVelocity(nbins, Rmax)
	radii_dm = [(_i+1)/float(nbins)*Rmax for _i in xrange(nbins)]


ax.set_xlabel('$x\ \mathrm{[kpc]}$')
ax.set_ylabel('$v\ \mathrm{[km\ s^{-1}]}$')

plt.subplots_adjust(left = 0.16, right = 0.98, bottom = 0.02, top = 0.98)

directory = "/home/rpverbek/articles/tbtf/sven2/"
fig.savefig(directory + 'compare_rcs.pdf'.format(sim, offset))


