import pytest
import numpy as np

from simulation.lambda_r import lambda_r, lambda_r_r_lim, lambda_r_profile, lambda_r_sim
from simulation.simdata import get_maps, get_phot

@pytest.fixture
def sim_name():
    return 'mb.71002_p100_a800_r600'

@pytest.fixture
def orbit_sideon():
    return True

@pytest.fixture
def maps(sim_name, orbit_sideon):
    maps = get_maps(sim_name, orbit_sideon=orbit_sideon)
    return maps

@pytest.fixture
def idx():
    return 35

@pytest.fixture
def one_map(maps, idx):
    return maps[idx]

@pytest.fixture
def bins():
    return 30

@pytest.fixture
def r_lim():
    return 5

@pytest.fixture
def center(sim_name, orbit_sideon, idx):
    phot = get_phot(sim_name, orbit_sideon=orbit_sideon).to_pandas()
    my_df = phot.loc[idx]
    center = (my_df.xcentroid, my_df.ycentroid)
    return center

@pytest.mark.array_compare
def test_lambda_profile(one_map, r_lim, bins):
    flux, v_los, sig_los = one_map['lum'], one_map['vlos'], one_map['sig_los']
    width = one_map.meta['WIDTH']
    return lambda_r_profile(flux, v_los, sig_los, width=width, r_lim=r_lim, bins=bins)

def test_lambda_r(one_map):
    flux, v_los, sig_los = one_map['lum'], one_map['vlos'], one_map['sig_los']
    assert lambda_r(flux, v_los, sig_los) == 0.18987114457712997

def test_lambda_r_center(one_map, center):
    flux, v_los, sig_los = one_map['lum'], one_map['vlos'], one_map['sig_los']
    assert lambda_r(flux, v_los, sig_los, center=center) == 0.19390877198811343

def test_lambda_r_r_lim(one_map, r_lim):
    flux, v_los, sig_los = one_map['lum'], one_map['vlos'], one_map['sig_los']
    width = one_map.meta['WIDTH']
    assert lambda_r_r_lim(flux, v_los, sig_los, width=width, r_lim=r_lim) == 0.18440312231219064
    # Test what happens if I include the whole image as with simple lambda_r()
    assert lambda_r_r_lim(flux, v_los, sig_los, width=width, r_lim=width*2) == lambda_r(flux, v_los, sig_los)

def test_lambda_r_r_lim_center(one_map, r_lim, center):
    flux, v_los, sig_los = one_map['lum'], one_map['vlos'], one_map['sig_los']
    width = one_map.meta['WIDTH']
    assert lambda_r_r_lim(flux, v_los, sig_los, width=width, r_lim=r_lim, center=center) == 0.1891922557707294
    assert lambda_r_r_lim(flux, v_los, sig_los, width=width, r_lim=r_lim, center=None) == 0.18440312231219064
    assert lambda_r_r_lim(flux, v_los, sig_los, width=width, r_lim=width*2, center=center) == lambda_r(flux, v_los, sig_los, center=center)


@pytest.mark.array_compare
def test_lambda_r_sim(maps):
    flux, v_los, sig_los = maps['lum'], maps['vlos'], maps['sig_los']
    return lambda_r_sim(flux, v_los, sig_los)

@pytest.mark.array_compare
def test_lambda_r_profiles_center(one_map, r_lim, bins, center):
    flux, v_los, sig_los = one_map['lum'], one_map['vlos'], one_map['sig_los']
    width = one_map.meta['WIDTH']
    return lambda_r_profile(flux, v_los, sig_los, width=width, r_lim=r_lim, bins=bins, center=center)


def test_lambda_profile_warning(one_map, r_lim):
    resolution = one_map.meta['RESOL']
    bins = resolution//2 + 1
    flux, v_los, sig_los = one_map['lum'], one_map['vlos'], one_map['sig_los']
    width = one_map.meta['WIDTH']
    with pytest.warns(None) as warning_list:
        lambda_r_profile(flux, v_los, sig_los, width=width, r_lim=r_lim, bins=bins)

    assert len(warning_list) == 2
    w = warning_list.pop()
    assert w.category == UserWarning
    assert str(w.message) == f'bins > resolution/2 ({bins} > {resolution//2}): this does not make a lot of sense...'

    # This is not important to be checked
    w = warning_list.pop()
    assert w.category == RuntimeWarning
    assert str(w.message) == f'invalid value encountered in double_scalars'
