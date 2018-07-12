import pynbody
import pynbody.filt as f
import numpy as np

path_to_snap = '/home/michele/sim/sim67001/snapshot_0036'
snap = pynbody.load(path_to_snap)

pynbody.analysis.halo.center(snap)
pynbody.analysis.angmom.faceon(snap.g, disk_size="2 kpc")
# r_eff_v = pynbody.analysis.luminosity.half_light_r(snap, band='v')
# print("half_light_r pynbody:", r_eff_v)

snap.s['feh'] = np.log10(snap.s['fest'] / snap.s['mass']) + 2.756433

# sub_snap = snap[snap['rxy'] < 1]
print(snap.s)
sub_snap = snap.s[np.logical_and.reduce((snap.s['feh'] > -5, snap.s['feh'] < 100, snap.s['rxy'] < 1))]
print(sub_snap)
Mstar = sub_snap['massform'].sum()
print("Mstar: ", Mstar.in_units(pynbody.units.Msol), "Msol")
