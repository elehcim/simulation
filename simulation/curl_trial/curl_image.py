from test_curl import M4_gradient, M4_gradient2D
import pynbody
from pynbody.plot.sph import image
import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':
    snap = pynbody.load('/home/michele/sim/Birmingham_dataset/snapshot_filaments/Visualization/m169p100s100.snap')
    pynbody.analysis.halo.center(snap.g)
    im = pynbody.sph.render_image(snap.g, 'rho', x2=4)
    plt.imshow(im)
    # image(snap.g, 'rho')
    # pynbody.sph.render_image(snap.g, 'vel', kernel=M4_gradient()

