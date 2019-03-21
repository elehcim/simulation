from astropy.modeling import models, fitting
from astropy.modeling.models import Sersic1D
import matplotlib.pyplot as plt
import matplotlib
import numpy as np


def order_of_magnitude(x):
    return np.floor(np.log10(x))

def compute_rms_error(x, y):
    return np.sqrt(np.sum((x-y)**2)/len(x))

def sersic_fit(sb_profile, r_eff, n_0, deviation_param=0.2, verblevel=0):
    r = sb_profile['rbins']
    sb = sb_profile['sb']
    idx = np.digitize(r_eff, r) - 1
    ampl = 10**(-0.4*sb)[idx]
    oom = 10**(order_of_magnitude(ampl))
    # print(oom)
    sersic_init = Sersic1D(amplitude=ampl/oom, r_eff=r_eff, n=n_0,
                           fixed={'amplitude': False,
                                  'r_eff': False,
                                  'n': False,
                                  },
                           # Allow wobbling a bit
                           bounds={'r_eff':(r_eff*(1-deviation_param), r_eff*(1+deviation_param)),
                                   #'amplitude':(ampl/oom*(1-deviation_param), ampl/oom*(1+deviation_param)),
                                   'n':(max(n_0-5, 0), n_0+5),  # That's the trick. I'm assuming something more or less continuos
                                  },
                          )

    # sersic_init.amplitude.min = 1e-11
    # sersic_init.n.min = 0.1
    # sersic_init.n.max = 10
    # print("Bounds: ", sersic_init.bounds)
    # fitter = fitting.LevMarLSQFitter()
    # fitter = fitting.SLSQPLSQFitter() # <--- this causes some error
    # The SimplexLSQFitter seems the best.
    fitter = fitting.SimplexLSQFitter()
    fitted_model = fitter(sersic_init, r, 10**(-0.4*sb)/oom, maxiter=1000)#, verblevel=verblevel)#, acc=1e-10)
    fitted_model.amplitude *= oom
    rms_error = compute_rms_error(sb, -2.5*np.log10(fitted_model(r)))
    # print(fitted_model)
    return fitted_model, fitter, rms_error

def plot_fitted_profile(p, fit):
    plt.figure(figsize=(8,5))
    plt.plot(p['rbins'], p['sb'], 'ko')
    plt.plot(p['rbins'], -2.5*np.log10(fit(p['rbins'])))
    plt.xlabel('R')
    plt.ylabel(r'$\mu_V$')
    plt.gca().invert_yaxis()
