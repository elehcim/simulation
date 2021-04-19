import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import argparse
from simulation.parsers.parse_trace import parse_trace
from simulation.units import gadget_time_units

def plot_adhoc(trace, norm=False, units='gadget'):
    if units == 'gyr':
        t = trace.t * gadget_time_units.in_units('Gyr')
        xlabel = 'Gyr'
    else:
        t = trace.t
        xlabel = ''
    if norm:
        plt.plot(t, trace.a_adhoc, label='$a_{\mathrm{adhoc}}$')
    else:
        plt.plot(t, trace.ax_adhoc, '--', label='$a_x$')
        plt.plot(t, trace.ay_adhoc, label='$a_y$')
        plt.plot(t, trace.az_adhoc, ':', label='$a_z$')
    # plt.plot(trace.t, trace.a, 'r', linewidth=1)
    plt.ylabel(r'$a_{\mathrm{adhoc}}/a_{tot}$')
    plt.xlabel(xlabel)
    plt.grid()
    plt.legend();
    plt.show()

def parse_args(cli=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(dest='sim_path', help="Path to the simulation snapshots")
    parser.add_argument('--units', choices=['gadget', 'gyr'], help="Units of x axis", default='gadget')
    parser.add_argument('--norm', help="Plot only norm", action='store_true')
    args = parser.parse_args(cli)
    return args


def main(cli=None):
    args = parse_args(cli)
    trace_path = os.path.join(args.sim_path, 'trace.txt')
    trace = parse_trace(trace_path)
    print(trace.head())
    plot_adhoc(trace, norm=args.norm, units=args.units)


if __name__=='__main__':
    main()
