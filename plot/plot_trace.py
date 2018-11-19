import numpy as np
import matplotlib.pyplot as plt
import sys

def main(filename=sys.argv[-1]):
    if not filename:
        print('cannot read input file')
        sys.exit(2)
    print('loading file...'.format())
    t = np.loadtxt(filename)
    print('done')
    plt.plot(t[:, 1], t[:,2])
    plt.plot((0,0),'r+')
    begin = t[0, 1], t[0, 2]
    end = t[-10, 1], t[-1, 2]
    plt.plot(*begin,'go')
    plt.plot(*end,'ro')
    plt.xlabel('x [kpc]')
    plt.ylabel('y [kpc]')
    plt.title('t: {:5.2f} - {:5.2f} Gyr'.format(t[0,0], t[-1,0]))
    plt.axes().set_aspect('equal')
    plt.show()

if __name__ == '__main__':
    main()