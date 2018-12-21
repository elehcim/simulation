import pynbody
from pynbody import filt
from .units import gadget_time_units, gadget_acc_units
import numpy as np

######################
## HYPLOT CODE  (hyplot/src/CParticle)
######################
# double CGasParticle::mgfe() const
# {
#   if (fe() != 0. && mg() != 0.)
#   	return log10(mg() / fe()) + 0.261299;
#   else
# 	return 0.471782;// -0.545184746714;
# }
# //______________________________________________________________________________
# /*!
#  * [Fe/H] = log10(Fe_mass / particleMass) + 2.756433 (see Valcke masters thesis).
#  * formula: [A/B] = log10(A/B)-log10(A_sol/B_sol)
#  *
#  * \return [Fe/H] for this particle. If Fe (or particle mass) is 0: -98.
#  */
# double CGasParticle::feh() const
# {
#   if (fe() != 0. && mass() != 0.)
#     return log10(fe() / mass()) + 2.756433;
#   else
#     return -98.;
# }






# snap.s['smooth'], snap.s['rho']

# snap.s['mgst_sph'] = snap.s.kdtree.sph_mean(snap.s['mgst'], 50)
# snap.s['fest_sph'] = snap.s.kdtree.sph_mean(snap.s['fest'], 50)
# snap.s['mass_sph'] = snap.s.kdtree.sph_mean(snap.s['mass'], 50)

# popIII_filt = filt.BandPass('feh', -5, 100)

MgFe_sol = -0.261299
FeH_sol = -2.756433

@pynbody.derived_array
def mgfe(snap):
    arr = np.log10(snap.s['mgst_sph']/snap.s['fest_sph']) - MgFe_sol
    arr[np.logical_or(snap.s['mgst_sph'] == 0.0, snap.s['fest_sph'] == 0.0)] = 0.471782
    return arr

@pynbody.derived_array
def feh(snap):
    arr = np.log10(snap.s['fest_sph']/snap.s['mass_sph']) - FeH_sol
    arr[np.logical_or(snap.s['fest_sph'] == 0.0, snap.s['mass_sph'] == 0.0)] = -98.0
    return arr



@pynbody.derived_array
def acce_norm(self):
    return np.sqrt((self['acce'] ** 2).sum(axis=1))

@pynbody.derived_array
def dt_acc(self, errtol=0.05, softening=0.03):
    return np.sqrt(2 * errtol * softening / self['acce_norm'])


#2 * All.CourantFac * SphP[p].Hsml / SphP[p].MaxSignalVel;
@pynbody.derived_array
def dt_courant(self, courant=0.1):
    return (2 * courant * self['smooth'] / self['cs']).in_units(gadget_time_units)

# @pynbody.derived_array
# def acce_norm(self):
#     arr = np.sqrt((self['acce'] ** 2).sum(axis=1))
#     arr.units = gadget_acc_units
#     return  arr

# @pynbody.derived_array
# def dt_acc(self, errtol=0.05, softening=0.03):
#     return (np.sqrt(2 * errtol * softening / self['acce_norm'] * pynbody.units.kpc)).in_units('kpc km**-1 s')