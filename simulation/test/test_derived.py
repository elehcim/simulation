import simulation
import pynbody
import pytest

@pytest.fixture
def snap():
    return pynbody.load('snapshot_0563_69p1')

@pytest.mark.array_compare
def test_feh(snap):
    arr = snap.g['feh']
    assert arr.units == pynbody.units.NoUnit()
    return arr

@pytest.mark.array_compare
def test_feh_stars(snap):
    arr = snap.s['feh']
    assert arr.units == pynbody.units.NoUnit()
    return arr

@pytest.mark.array_compare
def test_mgfe(snap):
    arr = snap.g['mgfe']
    assert arr.units == pynbody.units.NoUnit()
    return arr

@pytest.mark.array_compare
def test_mgfe_stars(snap):
    arr = snap.s['mgfe']
    assert arr.units == pynbody.units.NoUnit()
    return arr

@pytest.mark.array_compare
def test_neutral_fraction(snap):
    arr = snap.g['neutral_fraction']
    assert pynbody.units.Unit("1.00e+00")
    return arr


## TODO Don't forget stars