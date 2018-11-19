from __future__ import print_function
import argparse
import numpy as np
import matplotlib.pyplot as plt


from kick_ic import radial_period, plot_integrand, turnp2mom, NFWPotential, kpc_in_km



parser = argparse.ArgumentParser()

parser.add_argument('-m', '--halo-mass', dest='M_h', default=1e14, help='Mass of the analytical halo (Msol)')
parser.add_argument('-c', '--conc', dest='c', default=0, help='NFW halo concentration factor. Use 0 to compute c automatically using Wechsler2002 or Strigari2007 models', type=float)
parser.add_argument('-p', '--peri', dest='rp', default=200, help='Pericenter distance in kpc', type=float)
parser.add_argument('-a', '--apo', dest='ra', default=1000, help='Apocenter distance', type=float)
parser.add_argument('-r', '--radius', dest='r', default=800, help='Distance in kpc', type=float)

args = parser.parse_args()

rp = args.rp * kpc_in_km  # convert to km
ra = args.ra * kpc_in_km  # convert to km
r  = args.r  * kpc_in_km  # convert to km

nfw = NFWPotential(args.M_h, args.c)

E, J = turnp2mom(rp, ra, nfw.V0, nfw.dV0dr)
print("Binding energy = {:.2e} km^2/s^2".format(E))
print("J = {:.2e} km^2/s".format(J))

T_r, quad_error = radial_period(rp, ra, E, nfw.V0, J)

print("Radial period:    {} Gyr (err = {:.2g})".format(T_r, quad_error))
plot_integrand(rp, ra, E, nfw.V0, J)
plt.show()