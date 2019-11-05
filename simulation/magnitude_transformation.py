import numpy as np


# See jordi et al. 2006
# https://ui.adsabs.harvard.edu/abs/2006A%26A...460..339J/
# https://www.aanda.org/articles/aa/pdf/2006/46/aa6082-06.pdf


def color_sdss_g_r(VmR):
    """
    Get SDSS (g'-r') color from Johnson filters color (V-R)
    Jordi2006 eq. 7
    """
    return 1.646 * VmR - 0.139

def rmR(VmR):
    """
    Get (r'-R) color from Johnson filters color (V-R)
    Jordi2006 eq. 4
    """
    if VmR <= 0.93:
        return 0.267 * VmR + 0.088
    else:
        return 0.77 * VmR - 0.37

def get_sdss_r(V, R):
    VmR = V - R
    r = np.where(VmR <= 0.93,
        0.267 * VmR + 0.088,
        0.77 * VmR - 0.37) + R
    # r[VmR <= 0.93] = 0.267 * VmR + 0.088 + R
    # r[VmR > 0.93] = 0.77 * VmR - 0.37 + R
    return r

def color_sdss_r_i(RmI):
    """
    Get SDSS (r'-i') color from Johnson filters color (R-I)
    Jordi2006 eq. 2
    """
    return 1.007 * RmI - 0.236


def color_sdss_g_i(V, R, I):
    """
    Get SDSS (g'-i') color from Johnson filters V, R, I
    """
    return color_sdss_g_r(V-R) - color_sdss_r_i(R - I)