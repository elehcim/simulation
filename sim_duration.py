import os
import sys
from simulation import get_param_used
from snap_io import load_first_last_snap, snapshot_file_list
from datetime import datetime

# TODO these times (especially the timestamps) are UTC. We should convert it to local time

def unix2time(t):
    return datetime.fromtimestamp(t).strftime('%Y-%m-%d %H:%M:%S')

class SimDuration(object):
    """docstring for SimDuration"""
    def __init__(self, path):
        self.path = path
        f, l = load_first_last_snap(path)
        self.params = get_param_used(path)
        self.arrival = float(self.params['TimeMax'])
        self.f, self.l = f, l
        print('')
        self.f_creation = os.path.getmtime(f.filename)
        self.l_creation = os.path.getmtime(l.filename)
        self.dt = self.l_creation - self.f_creation
        self.dt_day = self.dt/3600/24
        self.gyr = l.properties['time'].in_units('Gyr') - f.properties['time'].in_units('Gyr')
        self._gyr_day = self.gyr/self.dt_day
        self._simtime_day = (self.l.header.time-self.f.header.time)/self.dt_day

    def __repr__(self):
        s =  "First created: {}\n".format(unix2time(self.f_creation))
        s += "Last created:  {}\n".format(unix2time(self.l_creation))
        s += "First time:    {:.4f} ({:.4f} Gyr)\n".format(self.f.header.time, self.f.properties['time'].in_units('Gyr'))
        s += "Last time:     {:.4f} ({:.4f} Gyr)\n".format(self.l.header.time, self.l.properties['time'].in_units('Gyr'))
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
        first = self.f if first is None else snaplist[first]
        last = self.l if last is None else snaplist[last]

        now = datetime.utcnow().timestamp()
        
        end = last.header.time
        begin = first.header.time

        eta = (self.arrival - begin) * self.dt / end
        return  unix2time(now + eta)



def gyr_day(path):
    f, l = load_first_last_snap(path)
    dt_day = (os.path.getmtime(l.filename) - os.path.getmtime(f.filename))/3600/24
    gyr = l.properties['time'].in_units('Gyr') - f.properties['time'].in_units('Gyr')
    return gyr/dt_day

