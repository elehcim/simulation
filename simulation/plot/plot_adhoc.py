import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import argparse
from simulation.parsers.parse_trace import parse_trace

def main(filename=sys.argv[-1]):
    arr = load_trace(filename)
    plot_trace_np(arr)
    plt.show()


def plot_adhoc(trace, trace_version):
    if trace_version == 1:
        ax_adhoc = trace.ax
        ay_adhoc = trace.ay
        az_adhoc = trace.az
    elif trace_version == 2:
        ax_adhoc = trace.ax_adhoc
        ay_adhoc = trace.ay_adhoc
        az_adhoc = trace.az_adhoc
    else:
        raise RuntimeError('Please specify a trace version 1 or 2')
    plt.plot(trace.t, ax_adhoc, '--',)
    plt.plot(trace.t, ay_adhoc)
    plt.plot(trace.t, az_adhoc, ':')
    # plt.plot(trace.t, trace.a, 'r', linewidth=1)
    plt.grid()
    plt.legend();
    plt.show()

def parse_args(cli=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(dest='sim_path', help="Path to the simulation snapshots")
    parser.add_argument('-v', '--trace-version', help="Version of the trace file", default=1, type=int)
    args = parser.parse_args(cli)
    return args


def main(cli=None):
    args = parse_args(cli)
    trace_path = os.path.join(args.sim_path, 'trace.txt')
    trace = parse_trace(trace_path, args.trace_version)
    print(trace.head())
    plot_adhoc(trace, args.trace_version)


if __name__=='__main__':
    main()
