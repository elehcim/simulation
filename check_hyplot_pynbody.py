import chyplot
import pynbody
import enums
import numpy as np
from numpy.testing import assert_array_equal, assert_allclose
import pytest

getProp = chyplot.cglobals.plmap.getSecond
"""
reader = chyplot.CDataGadget()
reader.setFilename('snapshot_0065')
data = reader.readFile()

s = pynbody.load('snapshot_0065')

# data.convertUnits() ## Actually it seems to do nothing??

# some common properties
x = np.array(data.getDataArray(enums.T_all, getProp('x'), True))
y = np.array(data.getDataArray(enums.T_all, getProp('y'), True))
z = np.array(data.getDataArray(enums.T_all, getProp('z'), True))
mass = np.array(data.getDataArray(enums.T_all, getProp('mass'), True))
v_chy = np.array(data.getDataArray(enums.T_all, getProp('velocity'), True))


birthtime = np.array(data.getDataArray(enums.T_star, getProp('birthtime'), True))
initialMass = np.array(data.getDataArray(enums.T_star, getProp('initialMass'), True))

# Fe = np.array(data.getDataArray(enums.T_star, getProp('[Fe/H]'), True))


# s.physical_units()  # This seems to cause problems because of the expansion factor
v_py_all = s['vel']

v_py = np.sqrt((v_py_all**2).sum(1))
np.testing.assert_array_equal(v_py, v_chy)

np.testing.assert_array_equal(mass, s['mass'])
np.testing.assert_array_equal(x, s['pos'][:,0])
np.testing.assert_array_equal(y, s['pos'][:,1])
np.testing.assert_array_equal(z, s['pos'][:,2])


np.testing.assert_array_equal(birthtime, s.s['tform'])
np.testing.assert_array_equal(initialMass, s.s['massform'])
# np.testing.assert_array_equal(Fe, s.g['Fe'])
"""
class Test_hyplot_pynbody(object):
	""" hyplot pynbody comparison"""

	def setup(self):
		reader = chyplot.CDataGadget()
		reader.setFilename('snapshot_0065')
		self.data = reader.readFile()
		# data.convertUnits() ## Actually it seems to do nothing??

		self.s = pynbody.load('snapshot_0065')

	# def setup_convert(self, test_mass_converted_units):
	# 	""" setup any state tied to the execution of the given method in a
	# 	class.  setup_method is invoked for every test method of a class.
	# 	"""
	# 	self.data.convertUnits()

	# def teardown_method(self, test_mass_converted_units):
	# 	""" teardown any state that was previously setup with a setup_method
	# 	call.
	# 	"""
	# 	self.data.unConvertUnits()

	def test_time(self):
		# These are the raw data from the Gadget snapshot ( s kpc km**-1)
		assert self.data.time() == float(self.s.properties['time'])

	def test_pos(self):
		x = np.array(self.data.getDataArray(enums.T_all, getProp('x'), True))
		y = np.array(self.data.getDataArray(enums.T_all, getProp('y'), True))
		z = np.array(self.data.getDataArray(enums.T_all, getProp('z'), True))
		assert_array_equal(x, self.s['pos'][:,0])
		assert_array_equal(y, self.s['pos'][:,1])
		assert_array_equal(z, self.s['pos'][:,2])

	def test_dist_to_z(self):
		assert_array_equal(np.array(self.data.getDataArray(enums.T_all, getProp('dist_to_z'), True)),
						   self.s['rxy'])

	def test_id(self):
		assert_array_equal(np.array(self.data.getDataArray(enums.T_all, getProp('id'), True)),
			self.s['iord'])

	def test_mass(self):
		mass = np.array(self.data.getDataArray(enums.T_all, getProp('mass'), True))
		assert_array_equal(mass, self.s['mass'].in_units("1e10 Msol"))
	
	# def test_mass_converted_units(self):
	# 	# self.data.convertUnits()
	# 	mass = np.array(self.data.getDataArray(enums.T_all, getProp('mass'), True))

	# 	assert_array_equal(mass, self.s['mass']/10.0**6) #.in_units("10**6 Msol"))
	# 	# self.data.unconvertUnits()

	def test_velocity(self):
		v_chy = np.array(self.data.getDataArray(enums.T_all, getProp('velocity'), True))
		v_py = np.sqrt((self.s['vel']**2).sum(1))
		assert_array_equal(v_py.in_units("km s**-1"), v_chy)

	def test_cylindrical_tangential_velocity(self):
		rot_vel = np.array(self.data.getDataArray(enums.T_all, getProp('rot_vel'), True))
		# diff = rot_vel - (-self.s['vcxy'])
		# max_diff = np.max(diff)
		# print(max_diff)  						# 0.0303299259394
		# max_diff_idx = np.where(diff == max_diff)
		# print(rot_vel[max_diff_idx]) 		    # 0.05414845
		# print(-self.s['vcxy'][max_diff_idx])  # 0.02381853
		assert_allclose(rot_vel, - self.s['vcxy'], rtol=1e-5, atol=3.1e-2)  # The sign is changed 

	def test_rho(self):
		assert_array_equal(np.array(self.data.getDataArray(enums.T_gas, getProp('density'), True)),
						   self.s.g['rho'].in_units("1e10 Msol kpc**-3"))

	def test_star_metallicity(self):
		Z = np.array(self.data.getDataArray(enums.T_star, getProp('metallicity'), True))
		assert_array_equal(Z, self.s.s['metals'])

	def test_internal_energy(self):
		assert_array_equal(np.array(self.data.getDataArray(enums.T_gas, getProp('spec_ener'), True)),
						   self.s.g['u'])

	def test_temperature(self):
		temp = np.array(self.data.getDataArray(enums.T_gas, getProp('temperature'), True))
		assert_array_equal(temp, self.s.g['temp'])

	def test_pressure(self):
		assert_array_equal(np.array(self.data.getDataArray(enums.T_gas, getProp('pressure'), True)),
						   self.s.g['pres'])

	def test_star_properties(self):
		assert_array_equal(np.array(self.data.getDataArray(enums.T_star, getProp('birthtime'), True)),
							self.s.s['tform'])
		assert_array_equal(np.array(self.data.getDataArray(enums.T_star, getProp('initialMass'), True)),
							self.s.s['massform'])

	@pytest.mark.skip(reason="In hyplot the current time is not passed to the particle visitor")
	def test_star_age(self):
		assert_array_equal(np.array(self.data.getDataArray(enums.T_star, getProp('age'), True)),
						   self.s.s['age'])

	def test_star_metals(self):
		assert_array_equal(np.array(self.data.getDataArray(enums.T_star, getProp('Fe'), True)),
						   self.s.s['fest'])
		assert_array_equal(np.array(self.data.getDataArray(enums.T_star, getProp('Mg'), True)),
						   self.s.s['mgst'])

	def test_gas_metals(self):
		assert_array_equal(np.array(self.data.getDataArray(enums.T_gas, getProp('Fe'), True)),
							self.s.g['fesp'])
		assert_array_equal(np.array(self.data.getDataArray(enums.T_gas, getProp('Mg'), True)),
							self.s.g['mgsp'])

	# def test_smoothing_length(self):
	# 	assert_array_equal(np.array(self.data.getDataArray(enums.T_gas, getProp('smoothingLength'), True)),
	# 						self.s.g['smooth'])
		# assert_array_equal(np.array(self.data.getDataArray(enums.T_star, getProp('smoothingLength'), True)),
		# 					self.s.s['smooth'])

class Test_hyplot_pynbody_unit_converted(object):
	""" hyplot pynbody comparison"""

	def setup(self):
		reader = chyplot.CDataGadget()
		reader.setFilename('snapshot_0065')
		self.data = reader.readFile()
		self.data.convertUnits()

		self.s = pynbody.load('snapshot_0065')

	def test_time(self):
		np.isclose(self.data.time(), self.s.properties['time'].in_units("Gyr"), rtol=1e-6)
#        		   8.14970841403679  8.149660006021044

	def test_pos(self):
		x = np.array(self.data.getDataArray(enums.T_all, getProp('x'), True))
		y = np.array(self.data.getDataArray(enums.T_all, getProp('y'), True))
		z = np.array(self.data.getDataArray(enums.T_all, getProp('z'), True))
		assert_allclose(x, self.s['pos'][:,0], rtol=1e-6)
		assert_allclose(y, self.s['pos'][:,1], rtol=1e-6)
		assert_allclose(z, self.s['pos'][:,2], rtol=1e-6)

	def test_dist_to_z(self):
		assert_allclose(np.array(self.data.getDataArray(enums.T_all, getProp('dist_to_z'), True)),
						self.s['rxy'], rtol=1e-6)

	def test_mass(self):
		mass = np.array(self.data.getDataArray(enums.T_all, getProp('mass'), True))
		assert_array_equal(mass, self.s['mass'].in_units("1e6 Msol"))
	
	def test_velocity(self):
		v_chy = np.array(self.data.getDataArray(enums.T_all, getProp('velocity'), True))
		v_py = np.sqrt((self.s['vel']**2).sum(1))
		assert_array_equal(v_py, v_chy)

	def test_cylindrical_tangential_velocity(self):
		rot_vel = np.array(self.data.getDataArray(enums.T_all, getProp('rot_vel'), True))
		# diff = rot_vel - (-self.s['vcxy'])
		# max_diff = np.max(diff)
		# print(max_diff)  						# 0.0303299259394
		# max_diff_idx = np.where(diff == max_diff)
		# print(rot_vel[max_diff_idx]) 		    # 0.05414845
		# print(-self.s['vcxy'][max_diff_idx])  # 0.02381853
		assert_allclose(rot_vel, - self.s['vcxy'], rtol=1e-5, atol=3.1e-2)

	def test_rho(self):
		# Need to multiply for a big factor otherwise comparing is troublesome
		assert_allclose(1e37*np.array(self.data.getDataArray(enums.T_gas, getProp('density'), True)),
						1e37*self.s.g['rho'].in_units("g cm**-3"), rtol=1e-4)

	def test_internal_energy(self):
		assert_array_equal(np.array(self.data.getDataArray(enums.T_gas, getProp('spec_ener'), True)),
						   self.s.g['u'].in_units('cm**2 s**-2'))

	def test_star_metallicity(self):
		Z = np.array(self.data.getDataArray(enums.T_star, getProp('metallicity'), True))
		assert_array_equal(Z, self.s.s['metals'])


	def test_temperature(self):
		temp = np.array(self.data.getDataArray(enums.T_gas, getProp('temperature'), True))
		assert_array_equal(temp, self.s.g['temp'])

	def test_pressure(self):
		assert_array_equal(np.array(self.data.getDataArray(enums.T_gas, getProp('pressure'), True)),
						   self.s.g['pres'])

	def test_star_birthtime(self):
		assert_allclose(np.array(self.data.getDataArray(enums.T_star, getProp('birthtime'), True)),
						   self.s.s['tform'].in_units("Gyr"), rtol=1e-6, atol=3e-5)

	def test_star_initialMass(self):
		assert_array_equal(np.array(self.data.getDataArray(enums.T_star, getProp('initialMass'), True)),
						   self.s.s['massform'].in_units("1e6 Msol"))

	@pytest.mark.skip(reason="Here the problem seems that we have different times")
	def test_star_age(self):
		assert_array_equal(np.array(self.data.getDataArray(enums.T_star, getProp('age'), True)),
						   self.s.s['age'])

