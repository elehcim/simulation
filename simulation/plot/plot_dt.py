import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from ..parsers.parse_info import parse_info

def plot_dt(df):
    """
    Plot timestep
    """
    plt.plot(df.step, df.dt)
    plt.xlabel('step')
    plt.ylabel('dt')
    plt.yscale('log', basey=2)

    ax_2 = plt.gca()
    ax_10 = ax_2.twinx()
    ax_time = ax_2.twiny()

    # print(ax1.get_ylim())
    ax_10.set_ylim(ax_2.get_ylim())
    ax_10.set_yscale('log', basey=10)

    ax_time.set_xlim(df.t.min(), df.t.max())
    plt.show()

def main(filename=sys.argv[-1]):
    print("Parsing file {} ...".format(filename))
    df = parse_info(filename)
    print("Plotting")
    plot_dt(df)


if __name__ == '__main__':
    main(sys.argv[-1])