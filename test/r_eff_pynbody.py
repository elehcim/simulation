import pynbody
import pynbody.filt as f
import numpy as np

path_to_snap = '/home/michele/sim/sim67001/snapshot_0036'
snap = pynbody.load(path_to_snap)

pynbody.analysis.halo.center(snap.s)

r_eff_v = pynbody.analysis.luminosity.half_light_r(snap, band='v')
print("half_light_r pynbody:", r_eff_v)


# sub_snap = snap[snap['r'] < r_eff_v]
# Mstar = sub_snap.star['mass'].sum()
# print("Mstar half_light_r function:", Mstar)

# sphere = f.Sphere(r_eff_v)
# sub_snap_sphere = snap.s[sphere(snap.s).view(np.bool)]
# Mstar_sphere = sub_snap_sphere.star['mass'].sum()
# print("Mstar sphere:", Mstar_sphere)
