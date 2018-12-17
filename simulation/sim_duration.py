import os
from simulation import get_param_used
from snap_io import load_first_last_snap, snapshot_file_list
from datetime import datetime
import pynbody
import argparse
import numpy as np

# TODO these times (especially the timestamps) are UTC. We should convert it to local time


def unix2time(t):
    """

    Parameters
    ----------
    t : int
        UNIX time

    Returns
    -------
    A string with the timestamp relative to the input UNIX time
    """
    return datetime.fromtimestamp(t).strftime('%Y-%m-%d %H:%M:%S')


class SnapTime(object):
    """Contain time values of a GadgetSnap"""
    def __init__(self, snap):
        """

        Parameters
        ----------
        snap : pynbody.SimSnap
            Snapshot
        """
        # self._snap = snap
        self.time = snap.header.time
        self.number = int(snap.filename[-4:])
        self.gyr = snap.properties['time'].in_units('Gyr')
        self.creation = os.path.getmtime(snap.filename)

    def __str__(self):
        return "{:7.4f} ({:7.4f} Gyr)".format(self.time, self.gyr)

    def __repr__(self):
        return "Snap time: {:.4f} ({:.4f} Gyr). Created {}".format(self.time, self.gyr, unix2time(self.creation))


class SimDuration(object):
    """Common utilities to compute duration of a Simulation"""
    def __init__(self, path, first=0, last=-1):
        """

        Parameters
        ----------
        path : str
            path of the simulation
        first : int
            first snapshot index
        last : int
            last snapshot index
        """
        self.path = path
        snaplist = snapshot_file_list(os.path.expanduser(self.path), include_dir=True)
        self.first_snap = pynbody.load(snaplist[first])
        self.last_snap = pynbody.load(snaplist[last])
        self.f = SnapTime(self.first_snap)
        self.l = SnapTime(self.last_snap)
        self.params = get_param_used(path)
        self.arrival = float(self.params['TimeMax'])
        self.gyr = self.l.gyr - self.f.gyr
        self.dt = self.l.creation - self.f.creation
        self.dt_day = self.dt/3600/24
        if self.dt_day == 0.0:
            self._gyr_day = np.inf
            self._simtime_day = np.inf
        else:
            self._gyr_day = self.gyr/self.dt_day
            self._simtime_day = (self.l.time-self.f.time)/self.dt_day

    def __repr__(self):
        s =  "First created: {}\n".format(unix2time(self.f.creation))
        s += "Last created:  {}\n".format(unix2time(self.l.creation))
        s += "First time:    {}  ({})\n".format(self.f, self.f.number)
        s += "Last time:     {}  ({})\n".format(self.l, self.l.number)
        s += "Arrival:       {}\n".format(self.arrival)
        s += "Gyr/day:       {:.4f} ({:.4f} Gyr)\n".format(self._simtime_day, self.gyr_day)
        s += "ETA:           {}\n".format(self.eta())
        try:
            s += "Particles (start): dm:{} g:{} s:{}\n".format( len(self.first_snap.dm), len(self.first_snap.g), len(self.first_snap.s) )
            s += "Particles (now):   dm:{} g:{} s:{}".format( len(self.last_snap.dm), len(self.last_snap.g), len(self.last_snap.s) )
        except AttributeError:
            pass
        return s

    @property
    def simtime_day(self):
        return self._simtime_day

    @property
    def gyr_day(self):
        return self._gyr_day

    def eta(self):
        """Estimated time of arrival"""
        f = self.f
        l = self.l

        eta = f.creation + (self.arrival - f.time) / (l.time - f.time) * (l.creation - f.creation)
        return  unix2time(eta)


def gyr_day(path):
    f, l = load_first_last_snap(path)
    dt_day = (os.path.getmtime(l.filename) - os.path.getmtime(f.filename))/3600/24
    gyr = l.properties['time'].in_units('Gyr') - f.properties['time'].in_units('Gyr')
    return gyr/dt_day

# def plot_sim_speed(path):
#     import matplotlib.pyplot as plt
#     import numpy as np
#     snaplist = snapshot_file_list(os.path.expanduser(path), include_dir=True)
#     times_map = map(os.path.getmtime, snaplist)
#     times = np.array(times_map)
#     plt.hist()


def main(cli=None):
    parser = argparse.ArgumentParser("Display some info on the simulation duration")
    parser.add_argument(help='Path of the snapshots', dest='path')
    parser.add_argument(default=0, type=int, nargs='?', help='First file (ordinal, default=0)', dest='first')
    parser.add_argument(default=-1, type=int, nargs='?', help='Last file (ordinal, default=-1)', dest='last')
    args = parser.parse_args(cli)
    print(SimDuration(args.path, args.first, args.last))


if __name__ == '__main__':
    main()
