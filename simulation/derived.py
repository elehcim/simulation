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
# //______________________________________________________________________________
# /*!
#  * log10(M_mgsol/M_fesol) = -0.261299, source: Grevesse et al. 2007 en 2010
#  * formula: [Mg/Fe] = log10(M_Mg/M_Fe)-log10(M_Mg_sol/M_Fe_sol)
#  *
#  * \return [Mg/Fe] for this particle
#  */
# double
# CStarParticle::mgfe() const
# {
#   return log10(mg() / fe()) + 0.261299;
# }
# //______________________________________________________________________________
# /*!
#  * [Fe/H] = log10(Fe_mass / particleMass) + 2.756433 (see Valcke masters thesis).
#  * formula: [A/B] = log10(A/B)-log10(A_sol/B_sol)
#  *
#  * \return [Fe/H] for this particle. If Fe (or particle mass) is 0: -98.
#  */
# double
# CStarParticle::feh() const
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

# To avoid log(0). np.nan are possible too, but checks when interpolating have to be changed too.
NA_VALUE_FEH = -98.0
NA_VALUE_MGFE = 0.471782

MgFe_corr = -0.261299  # from Hyplot: log10(M_mgsol/M_fesol) = -0.261299, source: Grevesse et al. 2007 en 2010
FeH_corr = -2.756433

def _get_ftype(snap):
    if snap._unifamily is pynbody.family.star:
        return 'st'
    elif snap._unifamily is pynbody.family.gas:
        return 'sp'
    else:
        return None


@pynbody.derived_array
def feh(snap, na_value=NA_VALUE_FEH):
    ftype = _get_ftype(snap)
    if not ftype:
        raise RuntimeError("Derived array 'feh' is available only for family gas or star")
    name = 'fe' + ftype
    arr = np.log10(snap[name]/snap['mass']) - FeH_corr
    arr[np.logical_or(snap[name] == 0.0, snap['mass'] == 0.0)] = na_value  # -98.0
    return arr


@pynbody.derived_array
def mgfe(snap, na_value=NA_VALUE_MGFE):
    ftype = _get_ftype(snap)
    if not ftype:
        raise RuntimeError("Derived array 'mgfe' is available only for family gas or star")
    name = 'mg' + ftype
    arr = np.log10(snap[name]/snap['fe' + ftype]) - MgFe_corr
    arr[np.logical_or(snap[name] == 0.0, snap['fe' + ftype] == 0.0)] = na_value # 0.471782
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

@pynbody.derived_array
def v_norm(self):
    return np.sqrt((self['vel'] ** 2).sum(axis=1))


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
def neutral_fraction(snap):
    if snap._unifamily is not pynbody.family.gas:
        raise RuntimeError("Derived array 'neutral_fraction' is available only for family gas")

    hi = get_HI_vec(snap['temp'],
            snap['feh'],
            snap['mgfe'],
            snap.ancestor.header.redshift,
            snap['rho'])
    return pynbody.array.SimArray(hi, units=pynbody.units.Unit("1"))


@pynbody.derived_array
def mass_HI(self):
    return self['neutral_fraction'] * self['mass']


@pynbody.derived_array
def rho_HI(self):
    return self['neutral_fraction'] * self['rho']


@pynbody.analysis.profile.Profile.profile_property
def sigma_HI(self):
    assert self.ndim is 2
    return (self['mass_HI']/self._binsize).in_units('Msol pc**-2')


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


# Magnitude derived array (conversion between Johnson and SDSS filters)

from .magnitude_transformation import get_sdss_u, get_sdss_g, get_sdss_r, get_sdss_i, get_sdss_z

def sdss_band_from_johnson_Jordi2006(s, band):
    if band == 'sdss_u':
        return get_sdss_u(s['u_mag'], s['b_mag'], s['v_mag'])
    if band == 'sdss_g':
        return get_sdss_g(s['b_mag'], s['v_mag'])
    if band == 'sdss_r':
        return get_sdss_r(s['v_mag'], s['r_mag'])
    if band == 'sdss_i':
        return get_sdss_i(s['r_mag'], s['i_mag'])
    if band == 'sdss_z':
        return get_sdss_z(s['r_mag'], s['i_mag'])
    else:
        raise NotImplementedError(f'SDSS band {band} has not been implemented yet')


def sdss_band_from_johnson_Jester2005(s, band):
    # http://www.sdss3.org/dr8/algorithms/sdssUBVRITransform.php
    if band=='sdss_r':
        r = s['v_mag'] + np.where(s['u_mag']-s['b_mag'] < 0,
                                  - 0.46 * (s['b_mag']-s['v_mag']) + 0.11,
                                  - 0.42 * (s['b_mag']-s['v_mag']) + 0.11)
        r = np.where(s['r_mag']-s['i_mag'] < 1.15, r, np.nan)
        return r
    elif band=='sdss_g':
        g = s['v_mag'] + np.where(s['u_mag']-s['b_mag'] < 0,
                                  0.64 * (s['b_mag']-s['v_mag']) - 0.13,
                                  0.60 * (s['b_mag']-s['v_mag']) - 0.12)
        g = np.where(s['r_mag']-s['i_mag'] < 1.15, r, np.nan)
        return g
    else:
        raise NotImplementedError(f'SDSS band {band} has not been implemented yet')

_CONVERSION_TABLES = {'Jordi2006'  : sdss_band_from_johnson_Jordi2006,
                      'Jester2005' : sdss_band_from_johnson_Jester2005}

_MY_TABLE = 'Jordi2006'

from pynbody.derived import lum_den_template
import functools

sdss_bands = ['sdss_u', 'sdss_g', 'sdss_r', 'sdss_i', 'sdss_z']

for band in sdss_bands:
    X = lambda s, b=str(band): _CONVERSION_TABLES[_MY_TABLE](s, band=b)
    X.__name__ = band + "_mag"
    X.__doc__ = band + " magnitude from analysis.luminosity.calc_mags and converted to SDSS filter using Jester 2005"""
    pynbody.derived_array(X)

    lum_den = functools.partial(lum_den_template,band)

    lum_den.__name__ = band + "_lum_den"
    lum_den.__doc__ = "Luminosity density in astronomy-friendly units: 10^(-0.4 %s_mag) per unit volume. " \
                      "" \
                      "The magnitude is taken from analysis.luminosity.calc_mags."%band
    pynbody.derived_array(lum_den)
