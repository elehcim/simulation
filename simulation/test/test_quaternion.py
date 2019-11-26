import numpy as np

import quaternion
# Rotation of -90deg around x axis
# The normalized is used in case I change the coordinates
quat = (1/np.sqrt(2) * np.quaternion(1, -1, 0, 0)).normalized()

# same rotation
# quat_bis = 1/np.sqrt(2) * np.quaternion(-1, 1, 0, 0)

vec = np.array([0,1,0])
rot1 = quaternion.rotate_vectors(quat, vec)
# rot2 = quaternion.rotate_vectors(quat_bis, vec)
# np.testing.assert_array_equal(rot1, rot2)

# Now, let's see if what I do is the same
v = np.hstack([np.zeros((len(vec[np.newaxis]), 1)), vec[np.newaxis]])
vq = quaternion.as_quat_array(v)
new_vec = quat * vq * quat.conj()
rot3 = np.squeeze(quaternion.as_float_array(new_vec)[:, 1:])

# Using directly the function
from derotate_simulation import rotate_vec
rot4 = np.squeeze(rotate_vec(vec[np.newaxis], quat))

# Testing if I'm doing the same with the function by hand or calling the function
np.testing.assert_allclose(rot3, rot4, err_msg="I'm not doing the same anymore")


# See if using quaternion.rotate_vectors is different from rotate_vec
np.testing.assert_allclose(rot1, rot3, err_msg="I have a problem")

# np.testing.assert_allclose(rot1, rot4, err_msg="This is weird")


# So we can conclude that doing quat * vq * quat.conj()
# is actually doing the rotation as in quaternion.rotate_vectors(quat, vec)
