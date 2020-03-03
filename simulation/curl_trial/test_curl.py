import math
import numpy as np
import threading
from pynbody.sph import Kernel, Kernel2D

class M4_gradient(Kernel):
    """Gradient of the cubic spline kernel (M4 in Price 2012)"""
    def __init__(self):
        super().__init__()

    def get_value(self, d, h=1):
        """Get the value of the kernel for a given smoothing length."""
        if d < 1:
            f = - 3. * d  + (9. / 4.) * d ** 2
        elif d < 2:
            f = - 0.75 * (2. - d) ** 2
        else:
            f = 0

        return f / (math.pi * h ** 3)


class M4_gradient2D(Kernel):
    def __init__(self, k_orig=M4_gradient()):
        self.h_power = 2
        self.max_d = k_orig.max_d
        self.k_orig = k_orig
        self.safe = threading.Lock()

    def get_value(self, d, h=1):
        import scipy.integrate as integrate
        import numpy as np
        return 2 * integrate.quad(lambda z: self.k_orig.get_value(np.sqrt(z ** 2 + d ** 2), h), 0, h)[0]


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    # r = np.linspace(0, 3, 100)
    q2 = np.arange(0, 4.01, 0.02)
    k = Kernel()
    g = M4_gradient()
    # NOTE the square root in the abscissa.
    plt.plot(q2**0.5, k.get_samples(), label=r'$W$')
    plt.plot(q2**0.5, g.get_samples(), label=r'$\nabla W$')
    # plt.scatter(q2**0.5, [k.get_value(x ** 0.5) for x in q2], marker='.', label=r'$W$')
    # plt.plot(q2**0.5, [g.get_value(x ** 0.5) for x in q2], label=r'$\nabla W$')
    plt.xlabel('q')
    plt.legend()

    plt.figure()
    plt.plot(q2**0.5)
    # k2 = Kernel2D()
    # g2 = M4_gradient2D()
    # plt.plot(q, k2.get_samples(), label=r'$W_{2D}$')
    # plt.plot(q, g2.get_samples(), label=r'$\nabla W_{2D}$')
    # plt.xlabel('q')
    # plt.legend()

    plt.show()