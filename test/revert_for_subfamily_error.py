import pynbody
s = pynbody.load("../snapshot_0065")
print(s)

print(s.g)

isinstance(s.g, pynbody.snapshot.SimSnap)
t1 = s.rotate_x(90)

t1.sim

t1._sim
t2=s.g.rotate_x(90)

print(t2.sim)

print(t2._sim)



t1.revert()
print("reverted t1")
t2.revert()