# From a Gadget run
# Begin Step 2586, Time: 9.85915, Systemstep: 0.000619746
# domain decomposition... 
# NTopleaves= 3844
# work-load balance=2.8704   memory-balance=2.04597
# exchange of 0000118666 particles
# domain decomposition done. 
# begin Peano-Hilbert order...
# Peano-Hilbert done.
# Start force computation...
# Tree construction.
# Tree construction done.
# Begin tree force.
# Debug using last particle
# P[-1].Pos  = (-7.0273e+02, -9.1753e+02, -1.0876e+02)
# nfw.center = (0.0000e+00, 0.0000e+00, 0.0000e+00)
# nfwd.acc_factor = 3.2512e+08
# nfwd.rho_s = 3.4707e-04
# nfwd.R_s   = 1.2012e+02
# GRAV_ACC[-1] = (-2.7715e+02, -2.8074e+02, 1.1207e+02)
# NFW_ACC[-1]  = (-2.1333e+02, -2.7854e+02, -3.3016e+01)
# Added NFW acceleration.
# tree is done.
# Start density computation...
# ngb iteration 1: need to repeat for 0000003762 particles.
# Start hydro-force computation...
# force computation done.
from __future__ import print_function, division

import numpy as np


Pos  = np.array((-7.0273e+02, -9.1753e+02, -1.0876e+02))

rho_s = 3.4707e-04
R_s   = 1.2012e+02
G = 43007.3
acc_factor = 4 * np.pi * G * rho_s * R_s**3
acc_factor_res = 3.2512e+08
print(acc_factor, acc_factor_res)
print(np.allclose(acc_factor, acc_factor_res, atol=1e-6))
r = np.linalg.norm(Pos)
fac = acc_factor/(r**3) * ( np.log(1 + r/R_s) - r / (R_s + r) );

out_acc = Pos * fac
print('NFW_ACC[-1]  = (-2.1333e+02, -2.7854e+02, -3.3016e+01')

print(out_acc)



Pos2 = np.array((-1.2562e+03, -1.3673e+03, -2.5250e+02))
out_acc = Pos2 * fac