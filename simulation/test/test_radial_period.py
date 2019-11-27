import simulation
from simulation.simdata import get_radial_period
import pynbody
import pytest

@pytest.fixture
def sim_name():
	return 'mb.69002_p200_a800_r600_new_adhoc'

def test_measured(sim_name):
	measured = get_radial_period(sim_name)
	computed = get_radial_period(sim_name, which='computed')
	assert measured == 3.2854537870035205
	assert computed == 3.2358367425354504
