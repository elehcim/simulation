import pynbody

gadget_time_units = pynbody.units.Unit('kpc km**-1 s')
gadget_vel_units = pynbody.units.Unit('km s**-1')
gadget_acc_units = gadget_vel_units/gadget_time_units
gadget_acc_units

amu = pynbody.units.NamedUnit("amu",1.660539e-27*pynbody.units.kg)

