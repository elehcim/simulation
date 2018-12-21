import pynbody
import numpy as np
from pynbody import units

# an IDL file is here:  but produces format 1 output


STD_PROPERTIES = dict(a=1.0,
                      z=0,
                      boxsize=11.111*units.kpc,
                      h=0.7,
                      omegaL0=0.72,
                      omegaM0=0.28,
                      time=0 * units.Unit("s kpc km**-1"))


def spherical_uniform_random_positions(radius):
    N = len(radius)
    phi = 2 * np.pi * np.random.sample(N)
    theta = np.arccos(np.random.sample(N)*2 - 1)
    r = radius.view(type=np.ndarray)
    x = r * np.cos(phi) * np.sin(theta)
    y = r * np.sin(phi) * np.sin(theta)
    z = r * np.cos(theta)
    return x, y, z


def box_uniform_random_positions(N, box_size, seed=None):
    np.random.seed(seed)
    x = box_size * np.random.sample(N)
    y = box_size * np.random.sample(N)
    z = box_size * np.random.sample(N)
    return x, y, z


# def get_pos(radius, dtype=np.float32):
#     gas_xyz = uniform_random_positions(radius)
#     return np.column_stack(gas_xyz).astype(dtype)


# def create_new_particles(N, box_size, tot_gas_mass, seed=None):
#     d = dict()
#     coords = box_uniform_random_positions(N, box_size, seed)
#     d['pos'] = coords
#     d['vel'] = np.zeros(coords.shape, dtype=f.g['vel'].dtype)
#     d['iord'] = np.arange(iord_max, iord_max + N) + 1
#     d['mass'] = np.ones(N, dtype=np.float32) * tot_gas_mass.to_value(10**10 * u.solMass) / N
#     d['fesp'] = np.zeros(N, dtype=np.float32)
#     d['mgsp'] = np.zeros(N, dtype=np.float32)
#     gadget_time_unit = u.kpc/u.km*u.s
#     # d['pres'] = profile.pressure.to_value(10**10 * u.solMass*u.km**2 /(u.kpc**3 * gadget_time_unit**2)).astype(np.float32)
#     d['rho'] = np.zeros(N, dtype=np.float32)  ## Seems that rho are not necessary (see clustep code (clustep.py:38))
#     d['smooth'] = np.ones(N, dtype=np.float32)
#     # d['temp'] = profile.temperature.value.astype(np.float32)
#     # d['u'] = profile.u.to_value((u.km/u.s)**2).astype(np.float32)
#     d['zsph'] = np.zeros(N, dtype=np.float32)
#     return d


def create_gas_box(N, box_size, tot_gas_mass, filename):
    s = pynbody.snapshot.new(gas=N)
    gas_xyz = box_uniform_random_positions(N, box_size)
    pos = np.column_stack(gas_xyz).astype(np.float32)
    s.g['pos'] = pynbody.array.SimArray(pos, units.kpc).astype(np.float32)
    s.g['vel'] = pynbody.array.SimArray(np.zeros_like(s.g['pos']), units.km/units.s).astype(np.float32)
    s.g['mass'] = pynbody.array.SimArray(np.ones(N) * tot_gas_mass / N, 10**10*units.Msol)
    print("You'll have a mass resolution of: {:.2g} Msol".format(tot_gas_mass/N))
    # in gCosmoOrb, the mass is not important since it is put to -1 and n Gadget
    # this is substituted with the typical (max particle mass of the snapshot) mass of the galaxy particle.

    # print(s.g['mass'])
    # gadget_time_unit = u.kpc/u.km*u.s
    # gadget_time_unit_pynbody = units.kpc/units.km*units.s
    # s.g['pres'] = pynbody.array.SimArray(profile.pressure.to_value(
    #               10**10 * u.solMass*u.km**2 /(u.kpc**3 * gadget_time_unit**2)),
    #               10**10 * units.Msol*units.km**2 /(units.kpc**3 * gadget_time_unit_pynbody**2 ).astype(np.float32)
    s.g['iord'] = pynbody.array.SimArray(np.arange(N)).astype(np.int32)
    zero_arr = pynbody.array.SimArray(np.zeros(N)).astype(np.float32)
    rest_of_fields = ['u', 'rho', 'smooth', 'fesp', 'mgsp', 'temp', 'zsph', 'p'] #  'tform'
    for name in rest_of_fields:
        s.g[name] = zero_arr

    # s.g['smooth'] = np.ones(N, dtype=np.float32)

    s.properties = STD_PROPERTIES
    s.properties['boxsize'] = box_size * pynbody.units.kpc
    print('writing ', filename)
    s.write(fmt=pynbody.snapshot.gadget.GadgetSnap, filename=filename)


# def get_all_keys(snap):
#     """return all the (non derived) keys for all the families"""
#     ak = set()
#     for fam in snap.families():
#         ak = ak.union(snap[fam].loadable_keys()).union(snap[fam].keys()).union(snap[fam].family_keys())
#     ak = [k for k in ak if not k in ["x", "y", "z", "vx", "vy", "vz"]]
#     ak.sort()
#     return ak

# This approach is not feasible because to relax the cluster should be evolved before being simulated
# together with the galaxy
# def attach_gas_particles(input_snap, cluster, profile):
#     f = input_snap
#     N = len(profile.radius)
#     s = pynbody.snapshot.new(dm=len(f.dm), gas=N+len(f.gas), star=len(f.star), order='gas,dm,star')
#     tot_gas_mass = cluster.gas.mass_enclosed(fc.dm.virial_radius)
#     new_part = create_new_particles(f, profile, tot_gas_mass)
#     for k in get_all_keys(f.gas):
#         print("Stacking ", k)
#         arr = f.g[k]
#         elements = new_part[k]
#         if k in ['pos', 'vel']:
#             stacked = pynbody.array.SimArray(np.vstack([arr, elements]), arr.units)
#         else:
#             stacked = pynbody.array.SimArray(np.hstack([arr, elements]), arr.units)
#         s.g[k] = stacked

#     for k in get_all_keys(f.dm):
#         s.dm[k] = f.dm[k]

#     for k in get_all_keys(f.s):
#         s.s[k] = f.s[k]

#     s.properties = f.properties.copy()
#     s.properties['z'] = f.properties['z']
#     # FIXME the resulting time is different!
#     return s


if __name__ == '__main__':
    N  = 18145
    tot_gas_mass = 2.77367036005671250e2  # Msol.  Same resolution as Loic's glass
    create_gas_box(N, 1.0, tot_gas_mass, 'pippo.glass')
