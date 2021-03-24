import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import argparse
from simulation.parsers.parse_trace import parse_trace

def plot_adhoc(trace):
    plt.plot(trace.t, trace.ax_adhoc, '--', label='$a_x$')
    plt.plot(trace.t, trace.ay_adhoc, label='$a_y$')
    plt.plot(trace.t, trace.az_adhoc, ':', label='$a_z$')
    # plt.plot(trace.t, trace.a, 'r', linewidth=1)
    plt.grid()
    plt.legend();
    plt.show()

def parse_args(cli=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(dest='sim_path', help="Path to the simulation snapshots")
    args = parser.parse_args(cli)
    return args


def main(cli=None):
    args = parse_args(cli)
    trace_path = os.path.join(args.sim_path, 'trace.txt')
    trace = parse_trace(trace_path)
    print(trace.head())
    plot_adhoc(trace)


if __name__=='__main__':
    main()
