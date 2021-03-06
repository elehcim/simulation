import os
import sys
import numpy as np
import matplotlib.pyplot as plt


def load_trace(filename):
    if not os.path.isfile(filename):
        print('cannot read input file')
        sys.exit(2)
    print('loading file...'.format())
    t = np.loadtxt(filename)
    print('done')
    return t


def plot_trace(x, y, t=None, ax=None):
    """
    Plot orbital path of the galaxy.
    Start is the green point, end is the red point.
    """
    if ax is None:
        fig, ax = plt.subplots()
    ax.plot(x, y)
    ax.plot((0,0),'r+')
    begin = x[0], y[0]
    end = x[-1], y[-1]
    ax.plot(*begin,'go', markersize=1)
    ax.plot(*end,'ro', markersize=1)
    ax.set_xlabel('x [kpc]')
    ax.set_ylabel('y [kpc]')
    if t is not None:
        ax.set_title('t: {:5.2f} - {:5.2f} Gyr'.format(t[0], t[-1]))
    ax.set_aspect('equal')
    # plt.show()


def plot_trace_np(arr):
    plot_trace(arr[:, 1], arr[:, 2], arr[:, 0])


def plot_trace_df(df, ax=None):
    plot_trace(df.x.values, df.y.values, df.t.values, ax=ax)


def _plot_trace_text_old(t):
    plt.plot(t[:, 1], t[:,2])
    plt.plot((0,0),'r+')
    begin = t[0, 1], t[0, 2]
    end = t[-1, 1], t[-1, 2]
    plt.plot(*begin,'go')
    plt.plot(*end,'ro')
    plt.xlabel('x [kpc]')
    plt.ylabel('y [kpc]')
    plt.title('t: {:5.2f} - {:5.2f} Gyr'.format(t[0,0], t[-1,0]))
    plt.axes().set_aspect('equal')
    plt.show()


def _plot_trace_df_old(df):
    plt.plot(df.x, df.y)
    plt.plot((0,0),'r+')
    begin = df.x.iloc[0], df.y.iloc[0]
    end = df.x.iloc[-1], df.y.iloc[-1]
    plt.plot(*begin,'go')
    plt.plot(*end,'ro')
    plt.xlabel('x [kpc]')
    plt.ylabel('y [kpc]')
    plt.title('t: {:5.2f} - {:5.2f} Gyr'.format(df.t.iloc[0], df.t.iloc[-1]))
    plt.axes().set_aspect('equal')
    plt.show()


def main(filename=sys.argv[-1]):
    arr = load_trace(filename)
    plot_trace_np(arr)
    plt.show()


if __name__ == '__main__':
    main(sys.argv[-1])