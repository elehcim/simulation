import simulation
import pynbody
import pytest

@pytest.fixture
def snap():
    return pynbody.load('snapshot_0563_69p1')

@pytest.mark.array_compare
def test_feh(snap):
    return snap.g['feh']

@pytest.mark.array_compare
def test_feh_stars(snap):
    return snap.s['feh']

@pytest.mark.array_compare
def test_mgfe(snap):
    return snap.g['mgfe']

@pytest.mark.array_compare
def test_mgfe_stars(snap):
    return snap.s['mgfe']

@pytest.mark.array_compare
def test_neutral_fraction(snap):
    return snap.g['neutral_fraction']


## TODO Don't forget stars