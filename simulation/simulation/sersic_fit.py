from astropy.modeling import models, fitting
from astropy.modeling.models import Sersic1D
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from astropy.modeling import Fittable1DModel
from astropy.modeling.parameters import Parameter
from collections import OrderedDict

def order_of_magnitude(x):
    return np.floor(np.log10(x))

def compute_rms_error(x, y):
    return np.sqrt(np.sum((x-y)**2)/len(x))

def sersic_fit_linear(sb_profile, r_eff, n_0, deviation_param=0.2, verblevel=0):
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


class MagnitudeSersic1D(Fittable1DModel):
    r"""
    One dimensional Sersic surface brightness profile in magnitude.

    Parameters
    ----------
    amplitude : float
        Surface brightness at r_eff.
    r_eff : float
        Effective (half-light) radius
    n : float
        Sersic Index.

    See Also
    --------
    Gaussian1D, Moffat1D, Lorentz1D

    Notes
    -----
    Model formula:

    .. math::

        \mu(r)=\mu_e + \frac{2.5 b_n}{\ln 10} \left[\left(\frac{r}{r_{e}}\right)^{(1/n)}-1\right]

    The constant :math:`b_n` is defined such that :math:`r_e` contains half the total
    luminosity, and can be solved for numerically.

    .. math::

        \Gamma(2n) = 2\gamma (b_n,2n)
    """
    mu_e = Parameter(default=1)
    r_eff = Parameter(default=1)
    n = Parameter(default=4)
    _gammaincinv = None

    @classmethod
    def evaluate(cls, r, mu_e, r_eff, n):
        """One dimensional Sersic profile function in magnitude scale."""
        if cls._gammaincinv is None:
            try:
                from scipy.special import gammaincinv
                cls._gammaincinv = gammaincinv
            except ValueError:
                raise ImportError('Sersic1D model requires scipy > 0.11.')

        k = 2.5/np.log(10) * cls._gammaincinv(2 * n, 0.5)
        return mu_e + k * ((r/r_eff) ** 1/n - 1)


    @property
    def input_units(self):
        if self.r_eff.unit is None:
            return None
        else:
            return {'x': self.r_eff.unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return OrderedDict([('r_eff', inputs_unit['x']),
                            ('mu_e', outputs_unit['y'])])


def sersic_fit_sb_profile_(sb_profile, **kwargs):
    r = sb_profile['rbins'].view(np.ndarray)
    sb = sb_profile['sb'].view(np.ndarray)
    return sersic_fit(r, sb, **kwargs)

def sersic_fit(r, sb, r_eff, n_0, deviation_relative_r_eff=0.2, deviation_abs_n=0.5, verblevel=0):
    idx = np.digitize(r_eff, r) - 1
    mu_e_init = sb[idx]
    sersic_init = MagnitudeSersic1D(mu_e=mu_e_init, r_eff=r_eff, n=n_0,
                                    fixed={'mu_e': False,
                                    'r_eff': False,
                                    'n': False,
                                    },
                                    # Allow wobbling a bit
                                    bounds={'r_eff':(r_eff*(1-deviation_relative_r_eff ), r_eff*(1+deviation_relative_r_eff)),
                                    #'mu_e':(mu_e_init*(1-deviation_param), mu_e_init*(1+deviation_param)),
                                    'n':(max(n_0-deviation_abs_n, 0), n_0+deviation_abs_n),  # That's the trick. I'm assuming something more or less continuos
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
    fitted_model = fitter(sersic_init, r, sb, maxiter=1000)#, verblevel=verblevel)#, acc=1e-10)
    rms_error = compute_rms_error(sb, fitted_model(r))
    # print(fitted_model)
    return fitted_model, fitter, rms_error


def plot_fitted_profile(p, fit):
    plt.figure(figsize=(8,5))
    plt.plot(p['rbins'], p['sb'], 'ko')
    plt.plot(p['rbins'], fit(p['rbins']))
    plt.xlabel('R [kpc]')
    plt.ylabel(r'$\mu_V$')
    plt.gca().invert_yaxis()
