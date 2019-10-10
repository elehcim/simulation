import warnings
from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np
import pynbody
from astropy import units as u



def lambda_r_formula(flux, r, v_los, v_disp):
    return np.sum(flux * r * np.abs(v_los)) / np.sum(flux * r * np.sqrt(v_los**2 + v_disp**2))


def lambda_r(flux, v_los, v_disp):
    resolution = flux.shape[-1]
    x = y = np.arange(-resolution/2, resolution/2)
    xx, yy = np.meshgrid(x, y)
    z = np.sqrt(xx**2 + yy**2)
    r = z
    return np.sum(flux * r * np.abs(v_los)) / np.sum(flux * r * np.sqrt(v_los**2 + v_disp**2))


def lambda_r_r_lim(flux, v_los, v_disp, width, r_lim):
    assert flux.ndim == v_los.ndim == v_disp.ndim == 2
    resolution = flux.shape[-1]
    # print(f"flux.shape={flux.shape}")
    x = y = np.arange(-resolution/2, resolution/2)
    xx, yy = np.meshgrid(x, y)
    z = np.sqrt(xx**2 + yy**2)
    r = z
    r_max = r_lim/width * resolution
    r[r>r_max] = 0
    return np.sum(flux * r * np.abs(v_los)) / np.sum(flux * r * np.sqrt(v_los**2 + v_disp**2))


def lambda_r_sim(flux, v_los, v_disp):
    assert flux.ndim == v_los.ndim == v_disp.ndim == 3
    resolution = flux.shape[-1]
    print(f"flux.shape={flux.shape}")
    x = y = np.arange(-resolution/2, resolution/2)
    xx, yy = np.meshgrid(x, y)
    z = np.sqrt(xx**2 + yy**2)
    r = z[np.newaxis]
    print(f"r.shape={r.shape}")
    return np.sum(flux * r * np.abs(v_los), axis=(1,2)) / np.sum(flux * r * np.sqrt(v_los**2 + v_disp**2), axis=(1,2))

def lambda_r_profile(flux, v_los, v_disp, width, bins=50):
    resolution = flux.shape[-1]
    if bins>resolution/2:
        warnings.warn(f'bins > resolution/2 ({bins} > {resolution//2}): this does not make a lot of sense...')
    prof = list()
    r_arr = np.linspace(0, width/2, bins)
    for r in r_arr:
        prof.append(lambda_r_r_lim(flux, v_los, v_disp, width=width, r_lim=r))
    return np.array(prof)