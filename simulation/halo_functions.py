import numpy as np
from astropy import cosmology
import astropy.units as u

cosmo = cosmology.Planck15
RHO_C = cosmo.critical_density0.to('solMass/kpc**3')

@u.quantity_input
def halo_Wechsler2002(M: u.solMass):
    '''Return the NFW concentration factor following Wechsler et al 2002 fitting formula
    M is in solar mass'''
    c = 20 * (M/(1e11 * u.solMass))**-0.13
    return c

@u.quantity_input
def halo_Strigari2007(M: u.solMass):
    '''Return the NFW concentration factor following Strigari et al 2007 fitting formula
    M is in solar mass.

    This formula is more suited for dwarf Galaxies (M < 1e8 Msol)'''
    c = 33 * (M/(1e8 * u.solMass))**-0.06
    return c

def halo_scaled_density(c, rho_c=RHO_C, overdensity_factor=200.0):
    rho_s = overdensity_factor * c*c*c * rho_c / (3 * np.log(1+c) - c/(1+c))
    return rho_s  # Msol/km^3

def halo_scaled_radius(M, c, rho_c=RHO_C, overdensity_factor=200.0):
    R_s = ((M / (4.0/3.0 * np.pi * overdensity_factor * rho_c)) ** (1.0 / 3.0)) / c
    return R_s

def halo_scaled_radius_direct(M):
    # Gentile 2004, eq. 9
    R_s = 5.7 * (M/(1e11 * u.solMass)) ** 0.46
    return R_s # kpc

def halo_virial_radius(M, rho_c=RHO_C, overdensity_factor=200.0):
    R_vir = ((M / (4.0/3.0 * np.pi * overdensity_factor * rho_c)) ** (1.0 / 3.0))
    return R_vir

@u.quantity_input
def compute_nfw_c(M: u.solMass):
    if M < 1e8 * u.solMass:
        print("Initializing NFW using Strigari2007 formula")
        c = halo_Strigari2007(M)
    else:
        print("Initializing NFW using Wechsler2002 formula")
        c = halo_Wechsler2002(M)
    return c
