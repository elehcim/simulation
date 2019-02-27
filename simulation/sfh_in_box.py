import simulation
import pynbody
from pynbody import filt
from itertools import tee
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import warnings


def sfh_snap(snap, massform=True, trange=False, bins=100, **kwargs):
    '''
    star formation history from pynbody.plot.stars

    By default, sfh will use the formation mass of the star.  In tipsy, this will be
    taken from the starlog file. Set massform=False if you want the final (observed)
    star formation history

    Parameters:
    ----------
    trange : list, array, or tuple
        Specifies the time range. len(t_range) must be 2.
    bins : int
        number of bins to use for the SFH.
    massform : bool
        decides whether to use original star mass (massform) or final star mass.

    Returns
    -------
    sfhist : pynbody.array.SimArray (Msol/yr)
        The SFR in each time bin
    thebins : pynbody.array.SimArray (Gyr)
        The bins edges
    '''
    if 'nbins' in kwargs:
        bins = kwargs['nbins']

    if ((len(snap.g)>0) | (len(snap.d)>0)):
        simstars = snap.star
    else:
        simstars = snap

    if trange:
        assert len(trange) == 2
    else:
        trange = [simstars['tform'].in_units("Gyr").min(), simstars['tform'].in_units("Gyr").max()]

    binnorm = 1e-9 * bins / (trange[1] - trange[0])

    trangefilt = filt.And(filt.HighPass('tform', str(trange[0]) + ' Gyr'),
                          filt.LowPass('tform', str(trange[1]) + ' Gyr'))
    tforms = simstars[trangefilt]['tform'].in_units('Gyr')

    if massform:
        try:
            weight = simstars[trangefilt]['massform'].in_units('Msol') * binnorm
        except (KeyError, pynbody.units.UnitsException):
            warnings.warn("Could not load massform array -- falling back to current stellar masses", RuntimeWarning)
            weight = simstars[trangefilt]['mass'].in_units('Msol') * binnorm
    else:
        weight = simstars[trangefilt]['mass'].in_units('Msol') * binnorm

    sfhist, thebins = np.histogram(tforms, weights=weight, bins=bins, **kwargs)

    return pynbody.array.SimArray(sfhist, "Msol yr**-1"), pynbody.array.SimArray(thebins, "Gyr")


def _pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def sfh(sim):
    """
    Returns the Star Formation Rate

    It computes the total mass of stars born between two snapshots.
    The new stars are spotted using the particle IDs
    which are new w.r.t. the previous snapshot.

    Returns:
    --------
    dt: np.ndarray
        the delta-time of each snapshot w.r.t. the previous. (Gyr)
    sfr: np.ndarray
        the mass of star formed between snapshot n and n-1. (Msol)
    """
    ns = list()
    ns_idx = list()
    snaps = sim
    for (a0, a1) in _pairwise(snaps):
        s0, s1 = a0.s['iord'].view(np.ndarray), a1.s['iord'].view(np.ndarray)
        new_stars = np.setdiff1d(s1, s0)
        ns_idx.append(np.where(np.isin(s1, new_stars))[0])
        ns.append(new_stars)

    # `dts` contains right borders of the time bin
    mf = [0.0]
    dts = [0.0]
    for (idx, (a0, a1)) in zip(ns_idx, _pairwise(snaps)):
        mf.append(np.sum(a1.s['massform'][idx].in_units('Msol')).view(np.ndarray))
        dts.append(a1.header.time - a0.header.time)

    massformed = np.array(mf)
    dt = np.array(dts)

    t_conv_fac = (u.kpc/(u.km/u.s)).to(u.yr)

    # from https://stackoverflow.com/a/37977222
    # The problem is that in the first place I have 0.0/0.0.
    # I'd like the first element of sfr to be 0.0
    # sfr = (massformed) / (dt * t_conv_fac)  # Msol/yr
    denom = dt * t_conv_fac
    sfr = np.divide(massformed, denom, out=np.zeros_like(massformed), where=(denom != 0.0))

    return dt, sfr


def bin_sfh(times, sfr, bins=100):
    trange = times.max() - times.min()
    dt_orig = trange / len(sfr)
    dt_new = trange / bins
    bin_width_ratio = dt_orig / dt_new
    hist, binedges = np.histogram(times, bins=bins, weights=sfr*bin_width_ratio)
    return hist, binedges


def compute_binned_sfh(sim, bins=100):
    dt, sfr = sfh(sim)
    hist, binedges = bin_sfh(sim.times, sfr, bins=bins)
    return hist, binedges


def plot_hist_sfh(hist, binedges, ax=None, drawstyle='steps', **kwargs):
    if ax is None:
        fig, ax = plt.subplots()
    # from here: https://stackoverflow.com/a/18611135
    left, right = binedges[:-1], binedges[1:]
    X = np.array([left, right]).T.flatten()
    Y = np.array([hist, hist]).T.flatten()
    ax.plot(X,Y, drawstyle=drawstyle, **kwargs)

    max_sfr = np.max(hist)
    if ax.get_ylim()[1] < 1.2 * max_sfr:
        ax.set_ylim(0.0, 1.2 * max_sfr)
    ax.set_xlabel('Time [Gyr]')
    ax.set_ylabel('SFR [Msol/yr]')


def plot_binned_sfh(sim, bins=100, ax=None, drawstyle='steps', **kwargs):
    if ax is None:
        fig, ax = plt.subplots()

    hist, binedges = compute_binned_sfh(sim, bins)

    # from here: https://stackoverflow.com/a/18611135
    left, right = binedges[:-1], binedges[1:]
    X = np.array([left, right]).T.flatten()
    Y = np.array([hist, hist]).T.flatten()
    ax.plot(X,Y, drawstyle=drawstyle, **kwargs)

    max_sfr = np.max(hist)
    if ax.get_ylim()[1] < 1.2 * max_sfr:
        ax.set_ylim(0.0, 1.2 * max_sfr)
    ax.set_xlabel('Time [Gyr]')
    ax.set_ylabel('SFR [Msol/yr]')
    return hist, binedges

def plot_binned_sfh_old(sim, bins=100, ax=None, in_gyr=True):
    if ax is None:
        fig, ax = plt.subplots()
    dt, sfr = sfh(sim)
    sfh_t = sim[0].header.time + np.cumsum(dt)

    if in_gyr:
        from simulation.units import gadget_time_units, gadget_dens_units
        conv_fac = gadget_time_units.in_units('Gyr')
        sfh_t *= conv_fac
        ax.set_xlabel('Time [Gyr]')
    else:
        ax.set_xlabel('Time')
    mybins = np.linspace(sfh_t[0], sfh_t[-1], bins)
    # binwidth = sfh_t[1] - sfh_t[0]
    sfhist, thebins, patches = ax.hist(sfh_t, weights=sfr, bins=bins, histtype='step')

    ax.set_ylabel('SFR [Msol/yr]')
    return ax

def plot_sfh(sim, ax=None):
    dt, sfr = sfh(sim)

    if ax is None:
        fig, ax = plt.subplots()

    # pynbody.plot.stars.sfh(sim.snap_list[-1], subplot=ax, bins=len(sim), trange=sim.t_range)
    # ax1.grid()
    ax.plot(sim.times, sfr)
    ax.set_xlabel('Time [Gyr]')
    ax.set_ylabel('SFR [Msol/yr]')
    max_sfr = np.max(sfr)
    if ax.get_ylim()[1] < 1.2 * max_sfr:
        ax.set_ylim(0.0, 1.2 * max_sfr)
    return ax

if __name__ == '__main__':

    SIM_PATH = '/home/michele/sim/MySimulations/hi_osc/mb.69002_p200_a800_r600/out'

    sim = simulation.Simulation(SIM_PATH)

    ns = list()
    ds = list()
    ns_idx = list()
    ds_idx = list()
    new_stars = list()
    died_stars = list()
    for (a0, a1) in _pairwise(sim):
        s0, s1 = a0.s['iord'].view(np.ndarray), a1.s['iord'].view(np.ndarray)
        print(a0, a1)
        new_stars = np.setdiff1d(s1, s0)
        died_stars = np.setdiff1d(s0, s1)
        ns_idx.append(np.where(np.isin(s1, new_stars))[0])
        ds_idx.append(np.where(np.isin(s0, died_stars))[0])
        ns.append(new_stars)
        ds.append(died_stars)


    # left border of the time bin
    time = list()
    mf = list()
    for (idx, (a0, a1)) in zip(ns_idx, _pairwise(sim)):
        time.append(a0.header.time)
        mf.append(np.sum(a1.s['massform'][idx]).view(np.ndarray))


    times = np.array(time)
    massformed = np.array(mf)

    conv_fac = (u.kpc/(u.km/u.s)).to(u.Gyr)
    dt, sfr = sfh(SIM_PATH)

    fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(12,4))

    pynbody.plot.stars.sfh(sim.snap_list[-1], subplot=ax1, bins=len(sim), trange=sim.t_range)
    # ax1.grid()
    ax2.plot(sim.times, sfr)
    ax2.set_ylim(0, None)
    # ax2.grid()
    ax2.set_xlabel('time')
    ax2.set_ylabel('SFR [Msol/yr]')