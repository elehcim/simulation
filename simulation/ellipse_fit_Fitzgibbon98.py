"""This code is DEPRECATED, use scikit-image instead"""

import numpy as np
from numpy.linalg import eig, inv

class EllipseFit():
    # from http://nicky.vanforeest.com/misc/fitEllipse/fitEllipse.html
    # which actually implements Fitzgibbon98 algorithm
    # But remember: https://stackoverflow.com/a/48002645
    # I modify the code to follow Jason indication that there should be only one positive eigenvalue
    # There is also another paper they say improving on Fitzgibbon: https://github.com/bdhammel/least-squares-ellipse-fitting/blob/master/media/WSCG98.pdf
    def __init__(self, x, y):

        x = x[:,np.newaxis]
        y = y[:,np.newaxis]
        D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
        S = np.dot(D.T,D)
        C = np.zeros([6,6])
        C[0,2] = C[2,0] = 2; C[1,1] = -1
        E, V =  eig(np.dot(inv(S), C))
        # n = np.argmax(np.abs(E))  # original
        n = np.argmax(E)
        a = V[:,n]
        self._a = a
        res1, res2 = self._ellipse_axis_length()
        if res1>=res2:
            self.a = res1
            self.b = res2
        else:
            self.a = res2
            self.b = res1

    @property
    def center(self):
        a = self._a
        b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
        num = b*b-a*c
        x0=(c*d-b*f)/num
        y0=(a*f-b*d)/num
        return np.array([x0,y0])

    @property
    def theta( self ):
        a = self._a

        b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
        if b == 0:
            if a > c:
                return 0
            else:
                return np.pi/2
        else:
            if a > c:
                return np.arctan(2*b/(a-c))/2
            else:
                return np.pi/2 + np.arctan(2*b/(a-c))/2

    def _ellipse_axis_length( self ):
        a = self._a

        b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
        up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
        down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
        down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
        res1=np.sqrt(up/down1)
        res2=np.sqrt(up/down2)
        return np.array([res1, res2])



def ef2patch(e, edgecolor='red'):
    import matplotlib.patches as mpatches
    return mpatches.Ellipse(e.center, 2*e.a, 2*e.b, e.theta*180/np.pi, edgecolor=edgecolor, facecolor='none')


def ef2EllipseModel(e):
    import astropy
    return astropy.models.Ellipse2D(amplitude=1, x_0=e.center[0], y_0=e.center[1], a=e.a, b=e.b, theta=e.theta)