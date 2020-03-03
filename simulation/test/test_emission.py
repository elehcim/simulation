import simulation
import pynbody
import pytest
import numpy as np

# RTOL = 3e-1

@pytest.fixture
def snap():
    return pynbody.load('snapshot_0563_69p1')

@pytest.fixture
def snap1():
    return pynbody.load('snapshot_0001_69p1')

# I'm testing with data from hyplot
@pytest.mark.array_compare(atol=1e-57, rtol=1e-1)
def test_cii_snap(snap):
    arr = snap.g['cii']
    print(np.min(arr[np.nonzero(arr)]))
    print(np.max(arr[np.nonzero(arr)]))

    assert arr.units == pynbody.units.Unit('erg s**-1 cm**-3')
    return arr

@pytest.mark.array_compare(atol=1e-57, rtol=1e-2)
def test_cii_snap1(snap1):
    arr = snap1.g['cii'][:100]
    print(np.min(arr[np.nonzero(arr)]))
    print(np.max(arr[np.nonzero(arr)]))
    assert arr.units == pynbody.units.Unit('erg s**-1 cm**-3')
    return arr

# @pytest.mark.array_compare
# def test_halpha(snap):
    # arr = snap.s['halpha']
    # assert arr.units == pynbody.units.Unit('erg s**-1 cm**-3')
    # return arr

