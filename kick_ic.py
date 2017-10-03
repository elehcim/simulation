#!/usr/bin/env python
# Script to put a galaxy in a certain orbit around a analytical DM halo.
# Uses code from Sven's Spiral galaxy calculator

# usage: kick_ic.py [-h] [-o OUTPUT_FILE] [-m MTOT] [-c C] [-p RP] [-a RA]
#                   [-r R]  input_file

# Inputs:
# - Input file (a Gadget2 snapshot)
# - Mass of the attactor
# - Pericenter distance
# - Apocenter Distance
# - Current distance from the attractor
# Output:
# - A snapshot file with the particles kicked with the desired velocity
#   in order to be in orbit (2 body orbit) around the attractor

# from __future__ import print_function
import argparse
import enums
import os
import sys
# from astropy import cosmology
# from astropy import units as u
import numpy as np

DEBUG_PARTICLE = 4

parser = argparse.ArgumentParser()
# parser.add_argument("--sim", "--simulation", default=60003, type=int)
# parser.add_argument('--simdir', default='~/sim', help='Path to the simultaions directory')
# parser.add_argument('-n', '--snapshot', '--snap', type=int, default=0, help='Snapshot')

# parser.add_argument('--newsim', default=0, help='A new simulation will be created. This is the name of the new simulation.')

parser.add_argument('input_file', help='Input file')

parser.add_argument('-o', '--output_file', '--out', help='Output file')

parser.add_argument('-m', '--halo-mass', dest='M_h', default=1e14, help='Mass of the analytical halo (Msol)')
parser.add_argument('-c', '--conc', dest='c', help='NFW halo concentration factor', type=float)

parser.add_argument('-p', '--peri', dest='rp', default=200, help='Pericenter distance in kpc', type=float)
parser.add_argument('-a', '--apo', dest='ra', default=1000, help='Apocenter distance', type=float)
parser.add_argument('-r', '--radius', dest='r', default=800, help='Distance in kpc', type=float)

                    #Time of pericenter passage in Gyr, after the current snapshot')
# parser.add_argument('-e', '--ecc', help='Eccentricity')

args = parser.parse_args()

import chyplot

rp = args.rp #* u.kpc # 15 # kpc. Pericenter distance
ra = args.ra #* u.kpc
r = args.r   #* u.kpc
c = args.c
M_h = args.M_h #* u.solMass
# prograde = not args.retrograde

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

kpc_in_km = 3.08568025e16  # conversion of kpc to km
gyr_in_s = 3.1556926e16  # conversion of gigayears to sec
# G = 1.32749351440e21  # Gravitional constant in km^3/(s^2 * 10^10 Msol)
G = 1.32749351440e11  # Gravitional constant in km^3/(s^2 * Msol)
Msol = 1.98855e30  # solar mass in kg
# Mpc_in_km = 3.08567e19  # km
G_SIunits = 6.67408e-11  # m^3 kg^-1 s^-2
## Cosmology part
# cosmo = cosmology.Planck15
# rho_c = cosmo.critical_density0.to('solMass/km**3')

def critical_density(z, h=0.67, omega_m=0.24):
    '''return the critical density at z in MSol/km**3'''
    h /= (kpc_in_km * 1000) # in 1/s
    return 3 * (h * 100)**2  /(8 * np.pi * G) * (omega_m * (1+z)**3 + 1 - omega_m)

rho_c = critical_density(0)
print rho_c  # Msol/km**3
print("Critical density {:.4g} Msol/km^3".format(rho_c))
# print("Critical density {:.4g} g/cm^3".format(rho_c/( Msol * 10**3 * 10**15)))
print("Critical density {:.4g} kg/m^3".format(rho_c * Msol/10**9))


def halo_Wechsler2002(M):
    '''Return the NFW concentration factor following Wechsler et al 2002 fitting formula
    M is in solar mass'''
    c = 20 * (M/1e8)**-0.13
    return c

def halo_Strigari2007(M):
    '''Return the NFW concentration factor following Strigari et al 2007 fitting formula
    M is in solar mass.

    This formula is more suited for dwarf Galaxies (M < 1e8 Msol)'''
    c = 33 * (M/1e8)**-0.06
    return c

def halo_scaled_density(c, rho_c=rho_c, overdensity_factor=200.0):
    rho_s = overdensity_factor * c*c*c * rho_c / (3 * np.log(1+c) - c/(1+c))
    return rho_s  # Msol/km^3

def halo_scaled_radius(M, c, rho_c=rho_c, overdensity_factor=200.0):
    R_s = ((M / (4/3 * np.pi * overdensity_factor * rho_c)) ** (1. / 3)) / c
    return R_s  # km

# G = 6.674e-11 m^3/(kg s^2)
# Msol=1.9891e30 kg
# G m^3/(kg s^2) * 10^10 * Msol kg * (1e9 (km^3/m^3))
# 1.32752534e21

rp *= kpc_in_km  # convert to km
ra *= kpc_in_km  # convert to km
r *= kpc_in_km  # convert to km

# TODO: to be put in a NFW class maybe
rho_s = halo_scaled_density(M_h, c)   # Msol/km^3
R_s = halo_scaled_radius(M_h, c)      # km 
print "Halo parameters:"
print "  Concentration factor    c =", c
print "  Halo scaled density rho_s = {:.2e} kg/m^3 ({:.2e} Msol/km^3)".format(rho_s * Msol/10**9, rho_s)
print "  Halo scaled radius    R_s = {:.2f} kpc".format(R_s/kpc_in_km)
# def V0(r):
#     return - 4 * np.pi * G * rho_s * R_s**3 * np.log(1 + r/R_s) / r

def V0(r, M=M_h, R_s=R_s, c=c):
    """Return the negative potential of a NFW density distribution
    From Annelies Cloet-Osselaert PhD thesis, appendix B. 
    I removed the constant (1/R_s)/(1+x_b) related to the fact that 
    the integration has been limited and not extended to infinity.
    Output units are in km^2/s^2

    """
    return G * M * (np.log(1+r/R_s)/r )/(np.log(1+c) - c/(1+c))

# def dV0dr(r):
#     return 4 * np.pi * G * rho_s * R_s**3 * ((r + R_s) * np.log(1 + r/R_s) - r) / (r**2 *(r+R_s))

def dV0dr(r, M=M_h, R_s=R_s, c=c):
    """Return the partial derivative of the negative potential of a NFW density distribution
    From Annelies Cloet-Osselaert PhD thesis, appendix B"""
    return G * M * (-np.log(1+r/R_s)/r**2 + 1/(r*R_s * (1+r/R_s))) / (np.log(1+c) - c/(1+c)) 


def turnp2mom(rm, rp, V0, dV0dr):
    """Convert turning points (pericenter and apocenter) to 
    binding energy E and angular momentum J.
    J must have same sign as rm !!
    """
    if abs(rm) < rp:
        print("Common case, V(rp)={:.4g}, V0(rm)={:.4g}".format(V0(rp), V0(rm)))
        E = (rp * rp * V0(rp) - rm * rm * V0(abs(rm))) / (rp * rp - rm * rm)  # km^2/s^2
        J = rp * rm * np.sqrt(2 * (V0(rp) - V0(abs(rm))) / (rm * rm - rp * rp))  # km/s

    elif abs(rm) == rp and rp != 0.0:
        J = rp ** 1.5 * np.sqrt(-dV0dr(rp)) if rm > 0 else -rp ** 1.5 * np.sqrt(-dV0dr(rp))
        E = V0(rp) - J * J / (2 * rp * rp)

    elif rp == 0.0:
        E = V0(rp)
        J = 0.0

    else:
        raise ValueError("rp ({}) should be less or equal than rm ({}). Quitting.".format(rp, rm))
    return E, J


def get_velocity(rp, ra, r):
    """Return the radial and azimuthal speed of an object orbiting a (rp,ra) orbit
    when it is at radius r"""
    if not rp < r < ra:
        raise ValueError("The following should be true: rperi < r < rapo ({} < {} < {})".format(rp, r, ra))

    E, J = turnp2mom(rp, ra, V0, dV0dr)
    print "Binding energy =", E
    print "J =", J 

    # E = V0 - r^2
    v = np.sqrt(2 * (V0(r) - E))

    print "velocity modulus", v, "km/s"
    
    v_theta = J / r
    v_r = np.sqrt(v**2 - v_theta**2)

    return v_r, v_theta

def polar_to_cartesian(x, y, v_r, v_theta):
    r = np.sqrt(x**2 + y**2)
    vx = v_r * x / r - v_theta * y / r
    vy = np.sqrt(v_r**2 + v_theta**2 - vx**2)
    return vx, vy

def get_particle_velocities(data):
    part_vx = np.array(data.getDataArray(enums.T_all, chyplot.cglobals.plmap.getSecond('vx'), True))
    part_vy = np.array(data.getDataArray(enums.T_all, chyplot.cglobals.plmap.getSecond('vy'), True))
    part_vz = np.array(data.getDataArray(enums.T_all, chyplot.cglobals.plmap.getSecond('vz'), True))
    return part_vx, part_vy, part_vz

def print_average_particle_velocities(data):
    part_vx, part_vy, part_vz = get_particle_velocities(data)
    print "Mean velocity vx={:.2f} vy={:.2f} vz={:.2f}".format(part_vx.mean(), part_vy.mean(), part_vz.mean())

def velocity_dispersion(v):
    return np.sqrt((v**2 - v.mean()**2).mean())

def data_velocity_dispersion(data):
    part_vx, part_vy, part_vz = get_particle_velocities(data)
    part_vx_disp = velocity_dispersion(part_vx)
    part_vy_disp = velocity_dispersion(part_vy)
    part_vz_disp = velocity_dispersion(part_vz)

    return part_vx_disp, part_vy_disp, part_vz_disp

v_r, v_theta = get_velocity(rp, ra, r)
print v_r, v_theta, "km/s"

# The galaxy is on the y axis
x = 0
z = 0
y = np.sqrt(r**2 - x**2 - z**2)  # still in km

vx, vy = polar_to_cartesian(x, y, v_r, v_theta)
vz = 0
print "vx = ", vx

position = np.array([x, y, z]) / kpc_in_km # to kpc

print "We are using file {}".format(args.input_file)

apsis = list(np.round(np.array([rp, ra, r])/kpc_in_km))
if args.output_file is None:
    gic_file = "{}.kicked_p{}_a{}_r{}_c{}.gic".format(os.path.basename(args.input_file), *(apsis + [c]))
else:
    gic_file = "{}".format(args.output_file)

print "The ICs file to be created {}".format(gic_file)
print "Galaxy position: ({:.2f}, {:.2f}, {:.2f}) kpc".format(*position)
print "Galaxy velocity: ({:.2f}, {:.2f}, {:.2f}) km/s".format(vx, vy, vz)

reader = chyplot.CDataGadget()
writer = chyplot.CWriteGadget()

# reader.setPrefix(os.path.join(os.path.expanduser(args.simdir), "sim%04.d"%simulation))
# reader.setRunNumber(simulation)
# reader.set_file(args.snapshot)

reader.setFilename(args.input_file)

try:
    data = reader.readFile()
except chyplot.IOError as e:
    print
    print "****Error reading file****"
    print args.input_file
    # print "simulation number:", simulation
    print e
    print e.what()
    sys.exit(12)

print "time: {:.2f} Gyr".format(data.time())
print "got file ", args.input_file
data.rcom(True, enums.T_star, 0, 0, 0, True)
data.vcom(True, enums.T_star)

print_average_particle_velocities(data)

v_disp = data_velocity_dispersion(data)

print "Velocity dispersion sigma_x={:.2f} sigma_y={:.2f} sigma_z={:.2f}".format(*v_disp)

print "Velocity kick in z direction: vz =", vz, "(It should be zero!)"

# debug 
print "Print properties of particle {} for debug purposes:".format(DEBUG_PARTICLE)
print "  Speed before kick    ({:.2f}, {:.2f}, {:.2f})".format(*data.findParticle(DEBUG_PARTICLE).velocity()), "km/s"
print "  Position before kick ({:.2f}, {:.2f}, {:.2f})".format(*data.findParticle(DEBUG_PARTICLE).position()), "kpc"

data.translate(enums.T_all, *position)  # Move the center of the galaxy to the calculated position
print "Kicking particles..."
data.kick(enums.T_all, vx, vy, vz)      # Change the velocity of the galaxy

print_average_particle_velocities(data)
v_disp_after = data_velocity_dispersion(data)
print "Velocity dispersion sigma_x={:.2f} sigma_y={:.2f} sigma_z={:.2f}".format(*v_disp_after)
print "  Speed after kick     ({:.2f}, {:.2f}, {:.2f})".format(*data.findParticle(DEBUG_PARTICLE).velocity()), "km/s" 
print "  Position after kick  ({:.2f}, {:.2f}, {:.2f})".format(*data.findParticle(DEBUG_PARTICLE).position()), "kpc"

outpath = os.path.join(os.getcwd(), gic_file)

try:
    # output_dir =  os.path.join(args.outdir, "ICs")
    # os.makedirs = output_dir
    print "Going to write: ", outpath
    writer.writeFile(data, outpath, enums.T_all)
except chyplot.IOError as e:
    print
    print "****Error writing file****"
    print outpath
    # print "simulation number:", newsim
    print e
    print e.what()
    sys.exit(13)

print "timeNewSimulation: {:.2f} Gyr".format(data.time())

## This code generates a .gic file. Now to do:
## change parameter file
## make directories on node
## Copy .gic file and parameterfile to cluster. 
## run as: mpirun -np numb_of_cores ./Gadget2 parameterfiles/namesim.param 2 > namesim.output 2> namesim.err & 
##          (the 2>namesim.err isn't necessary: it helps to determine what went wrong when the code crashes. 
## don't forget the final 2 this is a restart flag
