import simulation
import pynbody
from itertools import tee
import astropy.units as u
import numpy as np


def _pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def sfh(sim_path):
    sim = simulation.Simulation(sim_path)

    ns = list()
    ns_idx = list()
    new_stars = list()
    snaps = sim
    for (a0, a1) in _pairwise(snaps):
        s0, s1 = a0.s['iord'].view(np.ndarray), a1.s['iord'].view(np.ndarray)
        new_stars = np.setdiff1d(s1, s0)
        ns_idx.append(np.where(np.isin(s1, new_stars))[0])
        ns.append(new_stars)

    # time contains left border of the time bin
    mf = list()
    dts = list()
    for (idx, (a0, a1)) in zip(ns_idx, _pairwise(snaps)):
        mf.append(np.sum(a1.s['massform'][idx].in_units('Msol')).view(np.ndarray))
        dts.append(a1.header.time - a0.header.time)

    mf.append(0)
    dts.append(np.inf)

    massformed = np.array(mf)
    dt = np.array(dts)

    t_conv_fac = (u.kpc/(u.km/u.s)).to(u.yr)

    sfr = (massformed) / (dt * t_conv_fac)
    # lets finish with a zero:
    sfr = np.append(sfr, 0)
    dt = np.append(dt, 0)

    return dt, sfr


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

    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(nrows=2, figsize=(12,4))

    pynbody.plot.stars.sfh(sim.snap_list[-1], subplot=ax1, bins=len(sim), trange=sim.t_range)
    # ax1.grid()
    ax2.plot(sim.times, sfr)
    ax2.set_ylim(0, None)
    # ax2.grid()
    ax2.set_xlabel('time')
    ax2.set_ylabel('SFR [Msol/yr]')