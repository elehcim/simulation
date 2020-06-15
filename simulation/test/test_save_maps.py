import os
import pytest
import numpy as np
import pynbody
from simulation.save_maps import Imaging
import simulation.derived

@pytest.fixture
def snap():
    s = pynbody.load(os.path.join(os.path.dirname(__file__), 'snapshot_0563_69p1'))
    pynbody.analysis.halo.center(s.s)
    return s

@pytest.fixture
def orbit_sideon():
    return True

@pytest.fixture
def width():
    return 20

@pytest.fixture
def resolution():
    return 50

@pytest.fixture
def imaging(snap, width, resolution):
    return Imaging(snap, width, resolution)

@pytest.mark.array_compare
def test_imaging_sb_mag_map(imaging):
    return imaging.sb_mag()

@pytest.mark.array_compare
@pytest.mark.parametrize("band", ['v', 'sdss_r'])
def test_imaging_sb_lum_map(imaging, band):
    return imaging.sb_lum(band=band)


@pytest.mark.array_compare
def test_imaging_vlos_map(imaging):
    return imaging.v_los_map()

@pytest.mark.array_compare
def test_imaging_sigma_map(imaging):
    return imaging.v_disp_los_map()
