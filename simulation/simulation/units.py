import pynbody

gadget_length_units = pynbody.units.Unit('kpc')
gadget_vel_units = pynbody.units.Unit('km s**-1')
gadget_mass_units = pynbody.units.Unit('1e10 Msol')
# gadget_time_units = pynbody.units.Unit('kpc km**-1 s')
gadget_time_units = gadget_length_units/gadget_vel_units

gadget_acc_units = gadget_vel_units/gadget_time_units
gadget_dens_units = gadget_mass_units/gadget_length_units**3
gadget_angmom_units = gadget_mass_units * gadget_vel_units * gadget_length_units

amu = pynbody.units.NamedUnit("amu",1.660539e-27*pynbody.units.kg)

gravity_G = pynbody.array.SimArray(6.674e-11, pynbody.units.Unit("N kg**-2 m**2"))
