import pynbody
import matplotlib.pyplot as plt
import numpy as np
from simulation.derived import feh, mgfe


def mdf(snap, nbins=100, range=None, **kwargs):
    """
    Metallicity Distribution Function

    Parameters:
    -----------
    nbins : int
        Number of bins.
    range : 2-tuple
        default (-5, 0.3).
    kwargs : dict
        Passed to np.histogram

    Returns:
    --------
    metpdf : np.ndarray
        The probability density function relative to metallicity distribution function.
    bins : np.ndarray
        The array used to bin the particles.

    """
    if range is None:
        range = (-5, 0.3)
    metpdf, bins = np.histogram(snap['feh'], weights=snap['mass'],
                                bins=nbins, density=True, range=range, **kwargs)

    return metpdf, bins


def plot_mdf(snap, filename=None, clear=True, range=None, ax=False, nbins=100, **kwargs):
    """
    Metallicity Distribution Function
    from pynbody.analysis.metals

    Parameters:
    -----------

    Returns:
    --------
    None

    **Usage:**

    >>> import pynbody.plot as pp
    >>> pp.mdf(s.s,linestyle='dashed',color='k')
    """
    if ax:
        plt = ax
    else:
        import matplotlib.pyplot as plt
        if clear:
            plt.clf()
        plt.xlabel('[Fe / H]')
        plt.ylabel('PDF')

    metpdf, bins = mdf(snap, nbins=nbins, range=range)

    midpoints = 0.5 * (bins[:-1] + bins[1:])

    plt.plot(midpoints, metpdf, **kwargs)
    if (filename):
        print("Saving " + filename)
        plt.savefig(filename)

