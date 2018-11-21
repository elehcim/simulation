from snap_io import snapshot_file_list
from datetime import datetime
import pynbody
import argparse
import os
import numpy as np
import simulation
import asciiplotlib as apl

# # FIXME: decide whether to use header or snapshots
# def main(filename=None):
#     if filename is None:
#         filename = sys.argv[1]
#     fd = open(filename, "rb")
#     (name, length) = read_block_head(fd)
#     header_block = fd.read(256)
#     header = _construct_gadget_header(header_block)
#     print(header)


class SimParticle(object):
    """Common utilities to compute duration of a Simulation"""
    def __init__(self, path, first=0, last=-1, width=50, height=15):
        self.path = path
        sim = simulation.Simulation(os.path.expanduser(self.path))
        self.sim = sim
            #         s += "Particles (start): dm:{} g:{} s:{}\n".format( len(self.first_snap.dm), len(self.first_snap.g), len(self.first_snap.s) )
            # s += "Particles (now):   dm:{} g:{} s:{}".format( len(self.last_snap.dm), len(self.last_snap.g), len(self.last_snap.s) )
        self.width = width
        self.height = height

        self.dm =  [len(s.dm) for s in sim]
        self.gas = [len(s.g) for s in sim]
        self.s = [len(s.s) for s in sim]

    def plot_graph(self):
        plot_graph(self.sim.times, self.dm, self.gas, self.s, width=self.width, height=self.height)

    def plot_multigraph(self):
        plot_multigraph(self.sim.times, self.dm, self.gas, self.s, width=self.width, height=self.height)


def plot_graph(times, dm , gas, stars, width=50, height=15):
    fig = apl.figure()
    fig.plot(times, dm, label="dm", width=width, height=height)
    fig.plot(times, gas, label="g", width=width, height=height)
    fig.plot(times, stars, label="s", width=width, height=height)
    fig.show()


def plot_multigraph(times, dm , gas, stars, width=50, height=15):
    fig = apl.SubplotGrid((1,3), padding=(0, 0))

    fig[0,0].plot(times, dm, label="dm", width=width, height=height) # , symbol='points pt "o"'
    fig[0,1].plot(times, gas, label="g", width=width, height=height) # , symbol='points pt "o"'
    fig[0,2].plot(times, stars, label="s", width=width, height=height) # , symbol='points pt "o"'
    fig.show()


def main(cli=None):
    parser = argparse.ArgumentParser("Display some info on the simulation particles")
    parser.add_argument(help='Path of the snapshots', dest='path')
    parser.add_argument(default=0, type=int, nargs='?', help='First file (ordinal, default=0)', dest='first')
    parser.add_argument(default=-1, type=int, nargs='?', help='Last file (ordinal, default=-1)', dest='last')
    args = parser.parse_args(cli)
    p = SimParticle(args.path, args.first, args.last)
    # p.plot_graph()
    p.plot_multigraph()


if __name__ == '__main__':
    main()
