import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from parsers.parse_info import parse_info

def plot_dt(df):
    """
    Plot timestep
    """
    plt.plot(df.step, df.dt)
    plt.xlabel('step')
    plt.ylabel('dt')
    plt.yscale('log', basey=2)

    ax1 = plt.gca()
    ax2 = ax1.twinx()
    print(ax1.get_ylim())
    ax2.set_ylim(ax1.get_ylim())
    ax2.set_yscale('log', basey=10)

    plt.show()

def main(filename=sys.argv[-1]):
    df = parse_info(filename)
    plot_dt(df)


if __name__ == '__main__':
    main(sys.argv[-1])