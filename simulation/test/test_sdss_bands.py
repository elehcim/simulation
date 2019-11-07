import simulation
import pynbody
import pytest

@pytest.fixture
def snap():
	return pynbody.load('snapshot_0563_69p1')

def test_sdss_r(snap):
	# Jordi2006
	assert snap.s['sdss_r_mag'][12] == -3.2521589265520654
	# Jester2005
	# assert snap.s['sdss_r_mag'][12] == -3.0939065719412753


def test_sdss_bands(snap):
	assert snap.s['sdss_u_mag'][12] == -0.35435015753379484
	assert snap.s['sdss_g_mag'][12] == -2.3523495402141785
	assert snap.s['sdss_r_mag'][12] == -3.2521589265520654
	assert snap.s['sdss_i_mag'][12] == -3.6430877921246028

