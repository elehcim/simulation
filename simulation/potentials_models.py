import astropy.units as u
import numpy as np
from numba import jit

from simulation.halo_functions import RHO_C, compute_nfw_c


@jit(nopython=True)
def nfw_mass_enclosed_fast(r, r_s, rho_s):
    """Return the mass eclosed in solMass"""
    return 4*np.pi*rho_s*r_s**3*(np.log(1+r/r_s) - r/(r_s+r))


class NFWModel:
    """docstring for NFWModel"""
    @u.quantity_input
    def __init__(self, M: u.solMass):
        self.M = M
        self.c = compute_nfw_c(M)
        self.rho_s = self.halo_scaled_density()
        self.r_s = self.halo_scaled_radius()
        self._rho_s = self.rho_s.to_value("solMass/kpc3")
        self._r_s = self.r_s.to_value("kpc")

    def __repr__(self):
        return "NFW Model (M={:.2g}, c={:.2f}, rho_s={:.2g}, r_s={:.2g}".format(self.M, self.c, self.rho_s, self.r_s)

    def halo_scaled_radius(self, rho_c=RHO_C, overdensity_factor=200.0):
        r_s = ((self.M / (4.0/3.0 * np.pi * overdensity_factor * rho_c)) ** (1.0 / 3.0)) / self.c
        return r_s

    def halo_scaled_density(self, rho_c=RHO_C, overdensity_factor=200.0):
        rho_s = overdensity_factor * self.c**3 * rho_c / (3 * np.log(1+self.c) - self.c/(1+self.c))
        return rho_s

    @u.quantity_input
    def mass_enclosed(self, r: u.kpc):
        return 4*np.pi*self.rho_s*self.r_s**3*(np.log(1+r/self.r_s) - r/(self.r_s+r))

    def mass_enclosed_kpc(self, r):
        return 4*np.pi*self.rho_s*self.r_s**3*(np.log(1+r/self.r_s.to_value(u.kpc)) - r/(self.r_s.to_value(u.kpc)+r))

    def mass_enclosed_fast_python(self, r):
        """assuming r in kpc"""
        return 4*np.pi*self._rho_s*self._r_s**3*(np.log(1+r/self._r_s) - r/(self._r_s+r))

    def mass_enclosed_fast(self, r):
        """assuming r in kpc"""
        return nfw_mass_enclosed_fast(r, self._r_s, self._rho_s)

    @property
    def virial_mass(self):
        """Mass inside c*r_s"""
        return 4*np.pi*self.rho_s*self.r_s**3*(np.log(1+self.c)-(self.c/(1+self.c)))

    @property
    def virial_radius(self):
        """c*r_s"""
        return self.c * self.r_s

    def rho(self, r):
        return self.rho_s/((r/self.r_s) * (1+ r/self.r_s)**2)


def mass_cluster_enclosed(r, cluster_model=NFWModel(1e14*u.solMass)):
    """Return the mass enclosed within r kpc"""
    return cluster_model.mass_enclosed_fast(r)


def tidal_radius(r, galaxy_mass, cluster_model, peri, apo=800):
    """
    Return tidal radius in kpc

    Following King 1962
    http://articles.adsabs.harvard.edu/pdf/1962AJ.....67..471K
    Formula (11)
    """
    ecc = (apo - peri) / (apo + peri)
    mass_cl = cluster_model.mass_enclosed_fast(r)
    rt = r * (galaxy_mass / (mass_cl*(3 + ecc)))**(1.0/3.0)
    return rt
