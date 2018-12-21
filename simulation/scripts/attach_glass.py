import pynbody
import numpy as np

snap_name = "np_sim62002_snapshot_0067_p200_a800_r600_b60_rdw50.gic.ok"
glass_name = "ic18145.ic"
snap = pynbody.load(snap_name)

glass = pynbody.load(glass_name)

glass['mass'] = np.ones(len(glass['mass']), dtype=np.float32) * -1

f = pynbody.new(gas=len(snap.gas), star=len(snap.star), dm=len(snap.dm)+len(glass.gas))

# Adjust DM
for block in ['pos', 'vel', 'iord', 'mass']:
    f.dm[block] = np.concatenate((glass[block], snap.dm[block]))


for block, families in snap.get_block_list().items():
    print('block=', block)
    # if block in ['pos', 'vel', 'iord', 'mass']:
    #     f[block] = snap[block]
    for family in families:
        if family.name != 'dm':
            print('  family=', family.name)
            getattr(f, family.name)[block] = getattr(snap, family.name)[block]


f.properties = {**snap.properties, **dict(#time=snap.properties[np.float32(snap.header.time),
                                          boxsize=snap.header.BoxSize * pynbody.units.kpc,
                                          z=snap.header.redshift,
                                          a=1+1/snap.header.redshift,
                                          omegaL0=0.72,
                                          omegaM0=0.28,
                                          h=0.7)}
print(f.properties)
f.write(filename="snap_glass.ic", fmt=pynbody.gadget.GadgetSnap, use_time=True)
