import pynbody
s = pynbody.load("snapshot_0065")
s.properties['omegaM0'] = 0.28
redshift = pynbody.analysis.cosmology.redshift(s, s.properties['time'].in_units('Gyr'))
print("z (header) = {}".format(s.header.redshift))
print("computed z = {}".format(redshift))