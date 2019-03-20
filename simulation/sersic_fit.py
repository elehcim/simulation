from astropy.modeling import models, fitting
from astropy.modeling.models import Sersic1D
import matplotlib.pyplot as plt
import matplotlib
import numpy as np


def order_of_magnitude(x):
    return np.floor(np.log10(x))

def sersic_fit(profile, r_eff, n_0, deviation_param=0.1, verblevel=0):
    r = profile['rbins']
    sb = profile['sb']
    idx = np.digitize(r_eff, r) - 1
    ampl = 10**(-0.4*sb)[idx]
    oom = 10**(order_of_magnitude(ampl))
    # print(oom)
    sersic_init = models.Sersic1D(amplitude=ampl/oom, r_eff=r_eff, n=n_0,
                                  fixed={'amplitude': False,
                                         'r_eff': False,
                                         'n': False,
                                         },
                                  # Allow wobbling a bit
                                  bounds={'r_eff':(r_eff*(1-deviation_param), r_eff*(1+deviation_param)),
                                          #'amplitude':(ampl/oom*(1-deviation_param), ampl/oom*(1+deviation_param)),
                                         },
                                 )

    # sersic_init.amplitude.min = 1e-11
    sersic_init.n.min = 0.1
    sersic_init.n.max = 10
    # print("Bounds: ", sersic_init.bounds)
    fitter = fitting.SLSQPLSQFitter()
    s_fit = fitter(sersic_init, r, 10**(-0.4*sb)/oom, verblevel=verblevel)#, acc=1e-10)
    s_fit.amplitude *= oom
    # print(s_fit)
    return s_fit

def plot_fitted_profile(p, fit):
    plt.figure(figsize=(8,5))
    plt.plot(p['rbins'], p['sb'], 'ko')
    plt.plot(p['rbins'], -2.5*np.log10(fit(p['rbins'])))
    plt.xlabel('R')
    plt.ylabel(r'$\mu_V$')
    plt.gca().invert_yaxis()
