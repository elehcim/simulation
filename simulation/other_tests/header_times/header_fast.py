import numpy as np
import simulation
import timeit
import os
from time import time
from simulation import get_times_header, Simulation
from simulation.read_header import get_header
from simulation.snap_io import snapshot_file_list
SIM_DIR = "~/sim/MySimulations/ok_new_adhoc_or_not_affected/mb.69002_p200_a800_r600_new_adhoc/out"

def normal_way(sim_dir=SIM_DIR):
	sim = Simulation(sim_dir)
	# print(sim.times_header)
	return sim.times_header

def supposedly_faster_way(sim_dir=SIM_DIR):
	return get_times_header(sim_dir)

def get_times_header(sim_dir=SIM_DIR):
    """aster way to get times in the header"""
    snap_name_list = snapshot_file_list(os.path.expanduser(sim_dir), include_dir=True)
    return np.array([get_header(snap).time for snap in snap_name_list])

if __name__ == '__main__':

	start = time()
	# h1 = normal_way()
	s1 = time()
	# h2 = supposedly_faster_way()
	s2 = time()
	h3 = get_times_header()
	s3 = time()

	# print(f'supposedly_faster_way: {s2-s1} s')
	print(f'reading_header_only:   {s3-s2} s')
	# print(f'normal_way:            {s1-start} s')
	# np.testing.assert_equal(h2, h3)

	# np.testing.assert_equal(h1, h2)
	# print(h1, h2)
	number=10
	# timeit.timeit(normal_way, number=number)
	# timeit.timeit(supposedly_faster_way, number=number)