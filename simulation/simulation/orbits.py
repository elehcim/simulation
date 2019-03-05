import argparse
import os
import sys
# from astropy import cosmology
# from astropy import units as u
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
from scipy.integrate import simps
from simulation.units import gadget_time_units



kpc_in_km = 3.08568025e16  # conversion of kpc to km
gyr_in_s = 3.1556926e16  # conversion of gigayears to sec
# G = 1.32749351440e21  # Gravitional constant in km^3/(s^2 * 10^10 Msol)
G = 1.32749351440e11  # Gravitional constant in km^3/(s^2 * Msol)
Msol = 1.98855e30  # solar mass in kg
# Mpc_in_km = 3.08567e19  # km
G_SIunits = 6.67408e-11  # m^3 kg^-1 s^-2


def critical_density(z, h=0.67, omega_m=0.24):
    '''return the critical density at z in MSol/km**3'''
    h /= (kpc_in_km * 1000) # in 1/s
    return 3 * (h * 100)**2  /(8 * np.pi * G) * (omega_m * (1+z)**3 + 1 - omega_m)


RHO_C = critical_density(0)  # = 4.23932537713e-48 MSol/km**3

def halo_Wechsler2002(M):
    '''Return the NFW concentration factor following Wechsler et al 2002 fitting formula
    M is in solar mass'''
    c = 20 * (M/1e11)**-0.13
    return c


def halo_Strigari2007(M):
    '''Return the NFW concentration factor following Strigari et al 2007 fitting formula
    M is in solar mass.

    This formula is more suited for dwarf Galaxies (M < 1e8 Msol)'''
    c = 33 * (M/1e8)**-0.06
    return c


def halo_scaled_density(c, rho_c=RHO_C, overdensity_factor=200.0):
    rho_s = overdensity_factor * c*c*c * rho_c / (3 * np.log(1+c) - c/(1+c))
    return rho_s  # Msol/km^3


def halo_scaled_radius(M, c, rho_c=RHO_C, overdensity_factor=200.0):
    R_s = ((M / (4.0/3.0 * np.pi * overdensity_factor * rho_c)) ** (1.0 / 3.0)) / c
    return R_s  # km


def compute_c(M):
    if M < 1e8:
        c = halo_Strigari2007(M)
    else:
        c = halo_Wechsler2002(M)
    return c


class NFWPotential:

    def __init__(self, M_h, c=0.0):
        self.M = M_h
        if c == 0.0:
            print("Concentration factor 'c' computed automatically from halo mass ({:.1e} Msol)".format(M_h))
            self.c = compute_c(M_h)
        else:
            self.c = c
        self.R_s = halo_scaled_radius(self.M, self.c)    # km
        self.rho_s = halo_scaled_density(self.M, self.c) # Msol/km^3

    def V0(self, r):
        """Return the negative potential of a NFW density distribution
        From Annelies Cloet-Osselaert PhD thesis, appendix B.
        I removed the constant (1/R_s)/(1+x_b) related to the fact that
        the integration has been limited and not extended to infinity.
        Moreover constants in the potential do not affect almost anything.
        Output units are in km^2/s^2

        """
        return G * self.M * (np.log(1+r/self.R_s)/r)/(np.log(1+self.c) - self.c/(1+self.c))

# def dV0dr(r):
#     return 4 * np.pi * G * rho_s * R_s**3 * ((r + R_s) * np.log(1 + r/R_s) - r) / (r**2 *(r+R_s))

    def dV0dr(self, r):
        """Return the partial derivative of the negative potential of a NFW density distribution
        From Annelies Cloet-Osselaert PhD thesis, appendix B"""
        return G * self.M * (-np.log(1+r/self.R_s)/r**2 + 1/(r*self.R_s * (1+r/self.R_s))) / (np.log(1+self.c) - self.c/(1+self.c))


def turnp2mom(periapsis, apoapsis, V0, dV0dr):
    """Convert turning points (periapsis and apoapsis) to
    binding energy E and angular momentum J.

    V0 is the negative potential and dV0dr its derivative.

    Returns
    -------
    E in km^2/s^2
    J in km^2/s
    J must have same sign as rm !!
    """
    rm = periapsis
    rp = apoapsis
    if abs(rm) < rp:
        print("Common case, V(rp)={:.4g}, V0(rm)={:.4g}".format(V0(rp), V0(rm)))
        E = (rp * rp * V0(rp) - rm * rm * V0(abs(rm))) / (rp * rp - rm * rm)  # km^2/s^2
        J = rp * rm * np.sqrt(2 * (V0(rp) - V0(abs(rm))) / (rm * rm - rp * rp))  # km^2/s

    elif abs(rm) == rp and rp != 0.0:
        if dV0dr is None:
            raise ValueError("dV0dr should be provided in case rm==rp")
        J = rp ** 1.5 * np.sqrt(-dV0dr(rp)) if rm > 0 else -rp ** 1.5 * np.sqrt(-dV0dr(rp))
        E = V0(rp) - J * J / (2 * rp * rp)

    elif rp == 0.0:
        E = V0(rp)
        J = 0.0

    else:
        raise ValueError("rp ({}) should be greater than rm ({}). It's not. Quitting.".format(rp, rm))
    return E, J


def get_velocity(rp, ra, r, V0, dV0dr):
    """Return the radial and azimuthal speed of an object orbiting a (rp,ra) orbit
    when it is at radius r"""
    if not rp <= r <= ra:
        raise ValueError("The following should be true: rperi <= r <= rapo ({} <= {} <= {})".format(rp, r, ra))

    E, J = turnp2mom(rp, ra, V0, dV0dr)
    # print "Binding energy =", E, "km**2/s**2"
    # print "J =", J, "km**2/s"

    # E = V0 - v^2/2
    v = np.sqrt(2 * (V0(r) - E))

    print("Kick velocity modulus", v, "km/s")

    v_theta = J / r
    v_r = np.sqrt(v**2 - v_theta**2)

    return v_r, v_theta


def polar_to_cartesian(x, y, v_r, v_theta):
    r = np.sqrt(x**2 + y**2)
    vx = v_r * x / r - v_theta * y / r
    vy = np.sqrt(v_r**2 + v_theta**2 - vx**2)
    return vx, vy


def radial_period(r1, r2, E, V0, J):
    """Return the radial period.

    The time required for a particle in a spherically
    symmetric potential to travel from apocenter to pericenter and back.
    Source: GD sec: 3.1 eq. 3.17 (pag.146)

    Parameters
    ----------
    r1 : float
        Pericenter in km.
    r2 : float
        Apocenter in km.
    E : float
        Orbital energy in km**2/s**2.
    J : float
        Angular momentum in km**2/s.

    Returns
    -------
    T_r : float
        Radial period in Gyr.
    quad_err : float
        An estimate of the absolute error in the result.
    """

    if not r1 <= r2:
        raise ValueError("The lower integral limit r1 should be less than the upper one r2")

    def integrand(r):
        # Phi_r = -V0(r)
        return 2.0/gyr_in_s / np.sqrt(np.abs(2*(V0(r) - E) - J**2/r**2))

    # # With scipy.integrate.simps you can easily avoid the problematic points but the
    # # accuracy is much worse
    # x = np.linspace(rp, ra, 100000)[1:-1]
    # y = integrand(x)
    # res = simps(y,x)
    epsilon = 10
    T_r, quad_error = quad(integrand, r1+epsilon, r2-epsilon)  # Avoid apoasis because there 1/(dr/dt) is inf
    return T_r, quad_error


def radial_period_integrand(r, E, V0, J):
    """Binney and Tremaine 2008, Galactic Dynamics, sec: 3.1 eq. 3.17 (pag.146)

    Parameters
    ----------
    r : float
        Radius in km.
    E : float
        Orbital energy in km**2/s**2.
    V0 : function
        Spherially symmetric potential. Phi(r) = -V0(r). Output units are in km^2/s^2
    J : float
        Angular momentum in km**2/s.
    """
    res = 2.0 / np.sqrt(np.abs(2*(V0(r) - E) - J**2/r**2))
    return res


def plot_integrand(rp, ra, E, V0, J):
    # TODO fixme units
    fig = plt.figure()
    ax = fig.add_subplot(111)
    x_ = np.linspace(rp, ra, 100)[1:-1]
    print(len(x_))
    y_ = radial_period_integrand(x_, E, V0, J)
    ax.plot(x_,y_)
    # plt.show()


def measure_radial_period(sim):
    dr = np.diff(sim.r)
    zero_crossings = np.where(np.diff(np.signbit(dr)))[0]
    times_of_apsis = sim.times_header[zero_crossings]
    period_peri = times_of_apsis[2] - times_of_apsis[0]
    return period_peri * gadget_time_units.in_units('Gyr')


def compute_radial_period(rp, ra, pot=None):
    """Compute radial period around a spherical potential
    Parameters
    ----------
    rp : float
        Pericenter in kpc.
    ra : float
        Apocenter in kpc.
    pot : Potential
        Potential around which to compute the orbit. Default: NFWPotential(1e14, 0.0)

    Returns
    -------
    T_r : float
        Radial period in Gyr.
    quad_err : float
        An estimate of the absolute error in the result.
    """
    if nfw is None:
        M_h = 1e14
        pot = NFWPotential(M_h, 0.0)
    rp, ra, r = np.array(parse_simname(sim_name))*kpc_in_km
    E, J = turnp2mom(rp, ra, pot.V0, pot.dV0dr)
    T_r, quad_error = radial_period(rp, ra, E, pot.V0, J)
    print("Radial period:    {:.5} Gyr".format(T_r))
    return T_r, quad_error