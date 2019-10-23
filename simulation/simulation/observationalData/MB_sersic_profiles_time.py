
import matplotlib.pyplot as plt
import numpy as np
from hyplot.plot import PFigure

from matplotlib.ticker import ScalarFormatter
from hyplot.plot.PFormatDecorator import PSelectTickDecorator 
from hyplot.visual.PRunData import PRunData

from PyQt4.QtGui import QWidget

visual.fig = plt.figure(FigureClass = PFigure.PFigure)

ax = visual.fig.add_my_subplot(111)

#_______________________________________________________________________________
# plot the simulations 

sim     = 5015
labels  = [123]
colors  = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black', 'brown', 'purple']
for label, color in zip(labels, colors):
    data = PRunData(sim)
    x, profile, fit = data.getSBProfile(label)
    ax.plot(x, profile, '.', color=color, label=str(label*0.1))
    ax.plot(x, fit, color=color)

#_______________________________________________________________________________
# plot properties
ax.minorticks_on()

#ax.set_yscale('log')
#ax.set_xlim(-6, -24)
ax.set_ylim(35., 20)

fontdict = {'size':12}
# markerscale does not work : (
leg = ax.legend(loc="upper right", prop=fontdict, ncol=2, numpoints=1, shadow=True, fancybox=True, markerscale=2)

#for handle in leg.legendHandles:
#    print handle.get_xydata()

ax.set_xlabel('$r\ [\mathrm{kpc}]$')
ax.set_ylabel('$\mathrm{SB\ [mag/"}^2]$')

visual.finalize(name='DE_MB_sersic_profiles5015.png', show=False)
