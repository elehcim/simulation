# Script to put a galaxy in a certain orbit around a certain analytical halo.
# You can adjust it to put the position and velocity of the galaxy by hand.

import math
import os, sys, time
import argparse
import numpy as np
from scipy.optimize import newton

import enums

#Adjust!
gadgetDir = "/home/rpverbek/programs/gadget/"
trunkDir = "/home/rpverbek/programs/gadget/trunk/"


parser = argparse.ArgumentParser()
parser.add_argument("--sim", "--simulation", default=60003, type=int)
parser.add_argument('--newsim', default=0, help='A new simulation will be created. This is the name of the new simulation.')
parser.add_argument('--simdir', default='~/sim', help='Path to the simultaions directory')
parser.add_argument('--outdir', default='.', help='where to put the outputfile')


parser.add_argument('-n', '--snapshot', '--snap', type=int, default=0, help='Snapshot')
parser.add_argument('-m', '--halo-mass', dest='Mtot', default=1e12, help='Mass of the analytical halo')
parser.add_argument('-r', '--retrograde', default=False)

parser.add_argument('-p', '--rp', dest='r_p', default=15, help='Pericenter distance in kpc')
parser.add_argument('-t', '--tp', dest='t_p', default=2, help='Time of pericenter passage in Gyr, after the current snapshot')
parser.add_argument('-e', '--ecc', help='Eccentricity')


args = parser.parse_args()

import chyplot

r_p = args.r_p #15 # kpc. Pericenter distance
t_p = args.t_p #2 # Gyr. Time of pericenter passage, after the current snapshot
ecc = args.ecc
simulation = args.sim
prograde = not args.retrograde

# # **** Parameters to be set by hand ***
# # Simulation
# simulation = 14001
# # Snapshot
# snapshot = 36
# # Mass of the analytical halo
# Mtot = 1e12
# # This will create a new simulation. This is the name of the new simulation. 
# newsim = 1
# # Prograde or retrograde orbit
# prograde = True

# r_p = 15 # kpc. Pericenter distance
# t_p = 2 # Gyr. Time of pericenter passage, after the current snapshot
# ecc = 1  #Eccentricity. ecc<1 is an ellipse (ecc=0 is a circle), ecc=1 is a parabola, ecc>1 is a hyperbola.

kpc_in_km = 3.08568025e16 # conversion of kpc to km
gyr_in_s = 3.1556926e16 # conversion of gigayears to sec
G = 1.32749351440e21 #Gravitional constant in km^3/s^2/(10^10 Msol)

r_p *= kpc_in_km # convert to km
t_p *= gyr_in_s # convert to seconds

# Use this function in case of a parabola
def getPosVel_Parab(r_p, t_p, prograde, M_tot, printStuff=False):
	global G, kpc_in_km

	mu = M_tot*G # km^3/s^2 

	#Parameters of parabolic trajectory. See http://en.wikipedia.org/wiki/Parabolic_trajectory
	A = (3./2.)*math.sqrt(mu/(2*r_p**3))*t_p
	B = (A+math.sqrt(A**2+1))**1./3.
	nu = 2*math.atan(B-1/B)

	R = 2*r_p/(1+math.cos(nu))/kpc_in_km#r_p/math.cos(theta) / 3.08567758e16 Starting radius between 2 galaxies. Should preferably be sufficiently large

	theta = math.asin(math.sqrt(r_p/R/kpc_in_km))

	v = math.sqrt(mu*2./(R*kpc_in_km))
	v_p = math.sqrt(mu*2./(r_p))

	# Put printStuff to True if you want more output
	if printStuff:
		print 'nu =', nu/math.pi*180, 'degrees'
		print 'R =', R, 'kpc'
		print 'theta =', theta/math.pi*180, 'degrees'
		print 'v =', v, 'km/s'
		print 'v_p =', v_p, 'km/s'

	# x, y, z => distance where gas cloud should be placed (in kpc)
	x = -1*R*math.sin(nu)
	y = R*math.cos(nu)
	z = 0.

	# vx, vy, vz => velocities of the gas cloud (in km/sec = 1.0227 kpc/Gyr)
	vx = v*math.sin(math.pi-theta-nu)
	vy = v*math.cos(math.pi-theta-nu)
	vz = 0.

	# Change the x-position and -velocity in case of a retrograde orbit.
	if not prograde:
		x *= -1
		vx *= -1

	return x, y, z, vx, vy, vz

# Helper functions for the hyperbolic orbit
def get_time_to_peri(r, e, r_p, mu):
	#http://www.braeunig.us/space/orbmech.htm#hyperbolic
	a = r_p / (1-e)

	cosNu = (a * (1 - e**2) - r) / (e*r)

	coshF = (e + cosNu) / (1 + e*cosNu)
	F = np.arccosh(coshF)
	t = np.sqrt(-a**3 / mu) * (e * np.sinh(F) - F)

	return t

def f(r, t0, e, r_p, mu):
	return get_time_to_peri(r, e, r_p, mu) - t0

# Use this function in case of a hyperbola
def getPosVel_Hyperb(r_p, t_p, ecc, prograde, M_tot):
	global G

	mu = M_tot*G # km^3/s^2 

	Rstart = newton(f, 5*r_p, args=(t_p, ecc, r_p, mu))

	a = r_p / (1-ecc)
	cosNu = (a * (1 - ecc**2) - Rstart) / (ecc*Rstart)
	nu = np.arccos(cosNu)

	#https://en.wikipedia.org/wiki/Hyperbolic_trajectory

	x = Rstart*cosNu/kpc_in_km
	y = abs(Rstart*np.sin(nu))/kpc_in_km
	z = 0.

	if not prograde:
		y *= -1.

	v = np.sqrt( mu * (2. / Rstart - 1. / a) )

	if v > 1000:
		print "WARNING: very large velocity. v = {} km/s".format(v)

	#http://www.bogan.ca/orbits/kepler/orbteqtn.html

	phi = np.arctan( ecc * np.sin(nu) / (1 + ecc * cosNu) )

	vx = abs(v*np.cos(phi - nu + np.pi/2.))
	vy = abs(v*np.sin(phi - nu + np.pi/2.))
	vz = 0.

	if prograde:
		vy *= -1

	return x, y, z, vx, vy, vz

if ecc <= 1:
	if ecc < 1:
		print "Elliptical orbits aren't implemented yet..."
	print "Parabolic orbit"
	x, y, z, vx, vy, vz = getPosVel_Parab(r_p, t_p, prograde, args.Mtot)
else:
	print "Hyperbolic orbit"
	x, y, z, vx, vy, vz = getPosVel_Hyperb(r_p, t_p, ecc, prograde, args.Mtot)

#x, y, z, vx, vy, vz = 0, 0, 0, 0, 0, 0

# *************************************

# now, print the input stuff
print "We are using snapshot {} of sim {}.".format(args.snapshot, simulation)
gic_file = "sim{:05}.gic".format(args.newsim)
print "The ICs file to be created {}".format(gic_file)
print "The galaxy is put at distance: ({:.2f}, {:.2f}, {:.2f}) kpc and with speed of: ({:.2f}, {:.2f}, {:.2f}) km/s".format(x, y, z, vx, vy, vz)



totalData = None
reader = chyplot.CDataGadget()
writer = chyplot.CWriteGadget()

reader.setPrefix(os.path.join(os.path.expanduser(args.simdir), "sim%04.d"%simulation))
reader.setRunNumber(simulation)
reader.set_file(args.snapshot)

try:
	data = reader.readFile()
except chyplot.IOError as e:
        print
        print "****Error reading file****"
        print "simulation number:", simulation
        print e
        print e.what()
        sys.exit(12)

print "time: {:.2f} Gyr".format(data.time())
print "got sim ", simulation
data.rcom(True, enums.T_star, 0, 0, 0, True)
data.vcom(True, enums.T_star)

data.translate(enums.T_all,x, y, z) #Move the center of the galaxy to the calculated position
data.kick(enums.T_all, vx, vy, vz) #Change the velocity of the galaxy


try:
	# output_dir =  os.path.join(args.outdir, "ICs")
	# os.makedirs = output_dir
    writer.writeFile(data, os.path.join(args.outdir, gic_file), enums.T_all)
except chyplot.IOError as e:
    print
    print "****Error writing file****"
    print "simulation number:", newsim
    print e
    print e.what()
    sys.exit(13)

print "timeNewSimulation: {:.2f} Gyr".format(data.time())

## This code generates a .gic file. Now to do:
## change parameter file
## make directories on node
## Copy .gic file and parameterfile to cluster. 
## run as: mpirun -np numb_of_cores ./Gadget2 parameterfiles/namesim.param 2 > namesim.output 2> namesim.err & (the 2>namesim.err isn't necessary: it helps to determine what went wrong when the code crashes. 
## don't forget the final 2 this is a restart flag

