import pynbody
from pynbody import filt
from .units import gadget_time_units, gadget_acc_units
import numpy as np
import logging
import time

logger = logging.getLogger('simulation.derived')

######################
## HYPLOT CODE  (hyplot/src/CParticle)
######################
# double CGasParticle::mgfe() const
# {
#   if (fe() != 0. && mg() != 0.)
#       return log10(mg() / fe()) + 0.261299;
#   else
#   return 0.471782;// -0.545184746714;
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

MgFe_corr = -0.261299
FeH_corr = -2.756433

def _get_ftype(snap):
    if snap._unifamily is pynbody.family.star:
        return 'st'
    elif snap._unifamily is pynbody.family.gas:
        return 'sp'
    else:
        return None


@pynbody.derived_array
def feh(snap):
    ftype = _get_ftype(snap)
    if not ftype:
        raise RuntimeError("Derived array 'feh' is available only for family gas or star")
    name = 'fe' + ftype
    arr = np.log10(snap[name]/snap['mass']) - FeH_corr
    arr[np.logical_or(snap[name] == 0.0, snap['mass'] == 0.0)] = -98.0
    return arr


@pynbody.derived_array
def mgfe(snap):
    ftype = _get_ftype(snap)
    if not ftype:
        raise RuntimeError("Derived array 'mgfe' is available only for family gas or star")
    name = 'mg' + ftype
    arr = np.log10(snap[name]/snap['fe' + ftype]) - MgFe_corr
    arr[np.logical_or(snap[name] == 0.0, snap['fe' + ftype] == 0.0)] = 0.471782
    return arr


@pynbody.derived_array
def feh_smooth(snap):
    ftype = _get_ftype(snap)
    if not ftype:
        raise RuntimeError("Derived array 'feh' is available only for family gas or star")
    name = 'fe' + ftype
    name_sph = name + "_sph"
    nn = pynbody.config['sph']['smooth-particles']
    pynbody.sph.build_tree_or_trees(snap)
    snap.kdtree.set_array_ref('smooth',snap['smooth'])
    snap.kdtree.set_array_ref('mass',snap['mass'])
    snap.kdtree.set_array_ref('rho',snap['rho'])
    snap[name_sph] = snap.kdtree.sph_mean(snap[name], nn)
    snap['mass_sph'] = snap.kdtree.sph_mean(snap['mass'], nn)
    arr = np.log10(snap[name_sph]/snap['mass_sph']) - FeH_corr
    arr[np.logical_or(snap[name_sph] == 0.0, snap['mass_sph'] == 0.0)] = -98.0
    return arr


@pynbody.derived_array
def mgfe_smooth(snap):
    ftype = _get_ftype(snap)
    if not ftype:
        raise RuntimeError("Derived array 'mgfe' is available only for family gas or star")
    name = 'mg' + ftype
    name_sph = name + "_sph"
    nn = pynbody.config['sph']['smooth-particles']
    pynbody.sph.build_tree_or_trees(snap)
    snap.kdtree.set_array_ref('smooth',snap['smooth'])
    snap.kdtree.set_array_ref('mass',snap['mass'])
    snap.kdtree.set_array_ref('rho',snap['rho'])
    snap[name_sph] = snap.kdtree.sph_mean(snap[name], nn)
    snap['fest_sph'] = snap.kdtree.sph_mean(snap['fest'], nn)
    arr = np.log10(snap[name_sph]/snap['fest_sph']) - MgFe_corr
    arr[np.logical_or(snap[name_sph] == 0.0, snap['fest_sph'] == 0.0)] = 0.471782
    return arr



@pynbody.derived_array
def gas_metals(snap):
    # nn = pynbody.config['sph']['smooth-particles']
    # pynbody.sph.build_tree_or_trees(snap)
    # snap.kdtree.set_array_ref('smooth',snap['smooth'])
    # snap.kdtree.set_array_ref('mass',snap['mass'])
    # snap.kdtree.set_array_ref('rho',snap['rho'])
    # snap[name_sph] = snap.kdtree.sph_mean(snap[name], nn)
    # snap['mass_sph'] = snap.kdtree.sph_mean(snap['mass'], nn)
    arr = snap['fesp'] + snap['mgsp']/snap['mass']
    # arr[np.logical_or(snap[name_sph] == 0.0, snap['mass_sph'] == 0.0)] = -98.0
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

from .interp.gas_emission import get_HI_vec
@pynbody.derived_array
def HI(snap):
    if snap._unifamily is not pynbody.family.gas:
        raise RuntimeError("Derived array 'HI' is available only for family gas")

    hi = get_HI_vec(snap['temp'],
            snap['feh'],
            snap['mgfe'],
            snap.ancestor.header.redshift,
            snap['rho'])
    return hi


@pynbody.derived_array
def vz_disp(self):
    """SPH-smoothed local line-of-sigth velocity dispersion"""
    pynbody.sph.build_tree(self)
    nsmooth = pynbody.config['sph']['smooth-particles']
    self['rho']

    # logger.info(
    #     'Calculating velocity dispersion with %d nearest neighbours' % nsmooth)

    sm = pynbody.array.SimArray(np.empty(len(self['vz']),dtype=self['vz'].dtype),
                        self['vz'].units)

    self.kdtree.set_array_ref('rho',self['rho'])
    self.kdtree.set_array_ref('smooth',self['smooth'])
    self.kdtree.set_array_ref('mass',self['mass'])
    self.kdtree.set_array_ref('qty',self['vz'])
    self.kdtree.set_array_ref('qty_sm',sm)

    start = time.time()
    self.kdtree.populate('qty_disp',nsmooth)
    end = time.time()

    # logger.info('Velocity dispersion done in %5.3g s' % (end - start))

    return sm


@pynbody.derived_array
def mach(self):
    return (np.sqrt(self['v2']).in_units('km s**-1')/self['cs'].in_units('km s**-1'))
