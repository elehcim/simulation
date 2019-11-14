import numpy as np


# See jordi et al. 2006
# https://ui.adsabs.harvard.edu/abs/2006A%26A...460..339J/
# https://www.aanda.org/articles/aa/pdf/2006/46/aa6082-06.pdf


def get_sdss_u(U, B, V):
    g = get_sdss_g(B, V)
    u = g + 0.750 * (U-B) + 0.77*(B-V) + 0.72
    return u

def get_sdss_g(B, V):
    BmV = B - V
    g = 0.630 * BmV - 0.124 + V
    return g

def get_sdss_r(V, R):
    VmR = V - R
    r = R + np.where(VmR <= 0.93,
                     0.267 * VmR + 0.088,
                     0.77 * VmR - 0.37)
    return r

def get_sdss_i(R, I):
    i = I + 0.247 * (R-I) + 0.329
    return i

def get_sdss_z(R, I):
    z = R - 1.584 * (R-I) + 0.386
    return z

# Subaru SuprimeCam
def Aku_conversion_BmV_gmr_sec2_11_3(BmV):
    return (BmV - (0.2271+0.0038))/(1.3313-0.4216)
    # return (BmV - (0.2309))/(0.9097)

def Lupton2005_BmV_gmr(BmV):
    return (BmV - (0.2309))/(0.8914)


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
