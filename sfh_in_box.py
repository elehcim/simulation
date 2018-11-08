import simulation
from itertools import tee


def _pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def sfh(simpath):
    sim = simulation.Simulation(simpath)

    ns = list()
    ds = list()
    ns_idx = list()
    ds_idx = list()
    new_stars = list()
    died_stars = list()
    for (a0, a1) in _pairwise(sim):
        s0, s1 = a0.s['iord'].view(np.ndarray), a1.s['iord'].view(np.ndarray)
        new_stars = np.setdiff1d(s1, s0)
        died_stars = np.setdiff1d(s0, s1)
        ns_idx.append(np.where(np.isin(s1, new_stars))[0])
        ds_idx.append(np.where(np.isin(s0, died_stars))[0])
        ns.append(new_stars)
        ds.append(died_stars)


    # time contains left border of the time bin
    time = list()
    mf = list()
    for (idx, (a0, a1)) in zip(ns_idx, _pairwise(sim)):
        time.append(a0.header.time)
        mf.append(np.sum(a1.s['massform'][idx]).view(np.ndarray))


    times = np.array(time)
    massformed = np.array(mf)

    return times, massformed

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



