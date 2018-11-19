# trunk/scripts/render_gas.py
try:
    print "Importing mayavi..."
    from mayavi import mlab
    print "Done"
except:
    print "Please install mayavi to run this example"
    sys.exit(0)

import numpy as np
import chyplot
import enums

class Animation():
    def __init__(self, **kwargs):

        self.figsize = kwargs.get("figsize", (600,600))
        self.fig = mlab.figure(size=self.figsize)
        self.fig.scene.background = (1,1,1)
        self.fig.scene.foreground = (0,0,0)    

        self.reader = chyplot.CDataGadget()
        self.reader.setFilename("data/snapshot_0050")

        self.xlim = [-3, 3]
        self.ylim = self.xlim[:]
        self.zlim = self.xlim[:]

        self.renderPoints = kwargs.get("renderPoints", 50)

    def getData(self, filename=None, number=None):
        """
        Based on either a filename or a number: get the data.

        return : (s, sf, time)
        with s : gridded data
             sf: mlab scalar_field data
        """

        if filename:
            self.reader.setFilename(filename)
        else:
            self.reader.set_file(number)
        data = self.reader.readFile()
    
        # do rcom recentering
        data.rcom(True)
        data.convertUnits()
        
        gridder = chyplot.C3DSPHGridder(self.renderPoints,self.xlim[0], self.xlim[1], 
                                        self.ylim[0], self.ylim[1], self.zlim[0], self.zlim[1],
                                        enums.T_gas, data, "density")
        gridder.doGridding()
        s =  gridder.data()

        npt = self.renderPoints * 1j
        tx, ty, tz = np.mgrid[self.xlim[0]:self.xlim[1]:npt, self.ylim[0]:self.ylim[1]:npt, 
                              self.zlim[0]:self.zlim[1]:npt]
        sf = mlab.pipeline.scalar_field(tx, ty, tz, s)

        return s, sf, data.time()

    def plot(self, fname):
        s, sf, time = self.getData(fname)
        mlab.pipeline.iso_surface(
            sf, contours=[s.max()-0.9*s.ptp()], opacity=1.)
        
        #mlab.axes(nb_labels=5, line_width=2., xlabel="x", ylabel="y", zlabel="z",
        #          extent=(self.xlim[0], self.xlim[1], self.ylim[0], self.ylim[1], 
        #                  self.zlim[0], self.zlim[1]))

        mlab.savefig("figs/render_gas.png")
        mlab.show()

if __name__ == "__main__":
    anim = Animation()
    anim.plot("/home/michele/sim/MoRIA/sim71002/snapshot_0015")
    
