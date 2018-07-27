import os
import sys
from simulation import get_param_used
from snap_io import load_first_last_snap, snapshot_file_list
from datetime import datetime

# TODO these times (especially the timestamps) are UTC. We should convert it to local time

def unix2time(t):
    return datetime.fromtimestamp(t).strftime('%Y-%m-%d %H:%M:%S')


class SnapTime(object):
    """Contains time values of a GadgetSnap"""
    def __init__(self, snap):
        self._snap = snap
        self.time = snap.header.time
        self.number = int(snap.filename[:-4])
        self.gyr = snap.properties['time'].in_units('Gyr')
        self.creation = os.path.getmtime(snap.filename)

    def __str__(self):
        return "{:.4f} ({:.4f} Gyr)".format(self.time, self.gyr)

    def __repr__(self):
        return "Snap time: {:.4f} ({:.4f} Gyr). Created {}".format(self.time, self.gyr, unix2time(self.creation))


class SimDuration(object):
    """Common utilities to compute duration of a Simulation"""
    def __init__(self, path):
        self.path = path
        f, l = load_first_last_snap(path)
        self.params = get_param_used(path)
        self.arrival = float(self.params['TimeMax'])
        self.f, self.l = SnapTime(f), SnapTime(l)
        self.dt = self.l.creation - self.f.creation
        self.dt_day = self.dt/3600/24
        self.gyr = self.l.gyr - self.f.gyr
        self._gyr_day = self.gyr/self.dt_day
        self._simtime_day = (self.l.time-self.f.time)/self.dt_day

    def __repr__(self):
        s =  "First created: {}\n".format(unix2time(self.f.creation))
        s += "Last created:  {}\n".format(unix2time(self.l.creation))
        s += "First time:    {}  ({})\n".format(self.f, self.f.number)
        s += "Last time:     {}  ({})\n".format(self.l, self.l.number)
        s += "Arrival:       {}\n".format(self.arrival)
        s += "Gyr/day:       {:.4f} ({:.4f} Gyr)\n".format(self._simtime_day, self.gyr_day)
        s += "ETA:           {}".format(self.eta())
        return s

    @property
    def simtime_day(self):
        return self._simtime_day

    @property
    def gyr_day(self):
        return self._gyr_day

    def eta(self, first=None, last=None):
        """Estimated time of arrival"""
        snaplist = snapshot_file_list(os.path.expanduser(self.path), include_dir=True)
        first = self.f if first is None else SnapTime(snaplist[first])
        last = self.l if last is None else SnapTime(snaplist[last])

        now = datetime.utcnow().timestamp()

        end = last.time
        begin = first.time

        eta = (self.arrival - begin) * self.dt / end
        return  unix2time(now + eta)



def gyr_day(path):
    f, l = load_first_last_snap(path)
    dt_day = (os.path.getmtime(l.filename) - os.path.getmtime(f.filename))/3600/24
    gyr = l.properties['time'].in_units('Gyr') - f.properties['time'].in_units('Gyr')
    return gyr/dt_day

def main(arg=None):
    if arg is None:
        if len(sys.argv) < 2:
            raise ValueError("Please provide a simulation path")
        arg = sys.argv[1]
    print(SimDuration(arg))

if __name__ == '__main__':
    main()