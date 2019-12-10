import os
import glob
from simulation.util import get_aspect
import matplotlib.pyplot as plt
# from collection import namedtuple
import pprint

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data', 'AkuFigs')


def get_filename_from_y(y, data_dir=DATA_DIR):
    files = glob.glob(os.path.join(data_dir, '*.png'))
    # print(files)
    mapping = map(CMname, map(os.path.basename, files))
    # print(tuple(mapping))
    # d = dict()
    # for i, m in enumerate(mapping):
    #   d[m.y] = files[i]
        # print(m.y, files[i])
    d = dict((m.y, files[i]) for i, m in enumerate(mapping))
    # print(d)
    return d[y]

class CMname:
    def __init__(self, filename):
        self.filename = filename
        el = os.path.splitext(filename)[0].split('_')
        self.x = el[0]
        self.y = el[1]
        self.extent = tuple(map(float, el[2:]))

        self.elements = el

    def __repr__(self):
        return str(self.elements)

class AkuCM:
    def __init__(self, y):
        filename = get_filename_from_y(y)
        self.cmname = CMname(filename)
        self.extent = self.cmname.extent
        self.img = plt.imread(filename)
        self.aspect = get_aspect(self.extent, self.img.shape)

        self.filename = filename

    def __repr__(self):
        return pprint.pprint(self.__dict__)

if __name__ == '__main__':
    n = get_filename_from_y('gr')
    print(n)
    cm = AkuCM('mue')
    print(cm)