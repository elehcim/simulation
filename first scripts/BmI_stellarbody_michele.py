import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, LogLocator
import matplotlib.cm as cm

from hyplot.visual.PRunData import PRunData, PFileData
from hyplot.plot import PCanvas, PFigure, PPlot
from hyplot.gui  import PPlotDataOptions
from   hyplot import createColorBarFigs
import hyplot.comp.PGlobals         as glob
import chyplot
import math
import enums

from hyplot.plot import cplot, pythplot, PColorPlot, PScatterPlot, PCirclePlot, \
    quiver
from hyplot.comp import constants
import hyplot.comp.PGlobals as glob
import chyplot

from scipy.ndimage.filters import convolve, gaussian_filter
import numpy

matplotlib.rc('font', size=13)

solLumB = 5.441
solLumV = 4.820
solLumI = 4.08

globalsInitiated = False
def initGlobals():
    from hyplot     import PSettings
    from hyplot.gui import dataManager

    global globalsInitiated
    if not globalsInitiated:
        globalsInitiated = True
        
        glob.globals.settings    = PSettings.PSettings("HYPLOT", "Hyplot")
        glob.globals.dataManager = dataManager.PDataManager()
initGlobals()



### SETTINGS ###
SIM_PREFIX="/home/michele/sim/MoRIA/sim"
run = '62002'
sumrun = run
snap_first = 33
snap_last  = 45
offset = 0

gaussianFilter = True
sigma = 1

renderPlot = True

stellarPops = True
nrOfPops = 2      
constraints = ['age','age']#,'age']                             
limits = [[0.00,0.12],[0,0.01]]#,[0.05,0.1]]
colors = [(0.6,0.6,0),(1,1,0)]#,(0,0.6,0.6)]
size = [15,25]#,15]

# make plotargs dictionary

plotArgs = dict()

plotArgs['run']     = run
plotArgs['aspect']     = 'equal' 
plotArgs['xname']      = 'x'
plotArgs['yname']      = 'y'
plotArgs['xtitle']     = 'kpc'
plotArgs['ytitle']     = 'kpc'

plotArgs['xmin'], plotArgs['xmax'] = -3, 3
plotArgs['ymin'], plotArgs['ymax'] = -3, 3

plotArgs['plot1'] = dict()


# get analysis data from simulation

# sumdata = PRunData(sumrun)
# SF_data = sumdata.getData("SFR")
# time = sumdata.getData("time")


### PRODUCE PLOTS ###

print 'run',run
# glob.globals.dataManager.setMaximumData( 10 )

for snap in range(snap_first, snap_last+1):


    ### DO SETUP ###

    # set up hyplot canvas
    if renderPlot:
        canvas = PCanvas.create(setupMDI=False, showCanvas=False, width = 9, height = 5.5)
    else :
        canvas = PCanvas.create(setupMDI=False, showCanvas=False, **plotArgs)

    # read snapshot and set up data block
    dr = chyplot.CDataGadget(float(run))
    fdir = SIM_PREFIX+str(run)
    dr.setPrefix( fdir )
    dr.checkFilesPresent() # set the first and last dump
    dr.set_file( dr.firstDump() + snap)
    data = dr.readFile()
    print "time =", data.time()
    # pdata = glob.globals.dataManager.addData(data)

    # do stuff to blocks
    data.rcom(True, enums.T_star)
    data.vcom(True, enums.T_star)
    data.convertUnits()
    chyplot.cglobals.plmap.setDataBlock(data)
    
    # set up plotargs 
    plotArgs['fileNum'] = snap
    # plotArgs['filename'] = glob.globals.dataManager.activePData().filename
    # plotArgs['dataNumber'] = pdata.number

    # set up axes
    if renderPlot:
        axes   = canvas.figure.add_my_subplot(121)
        renderAxes = canvas.figure.add_my_subplot(122)
    else :
        axes   = canvas.figure.add_my_subplot(111)


    ### B-I COLOR PLOT ###

    plotArgs['cmin'], plotArgs['cmax'] = 0.5, 1.8001 #1.5
    plotArgs['plot1']['plotType']  = chyplot.CPyPlot.PLOT_GRID
    plotArgs['plot1']['nx']        = 150
    plotArgs['plot1']['ny']        = 150
    plotArgs['plot1']['type']      = 4
    plotArgs['plot1']['colormap']  = cm.get_cmap("RdBu_r")#("pink")
    plotArgs['plot1']['average']   = False
    plotArgs['plot1']['logscale']  = False
    plotArgs['plot1']['pointSize'] = 1

    # producing luminosity grids

    plotArgs['binname'] = 'luminosity_B'
    axes.setRanges(**plotArgs)
    axes.setProperties(**plotArgs)
    newcplotB = cplot.createCPlot(axes, **plotArgs['plot1'])  
    PyPlotB = pythplot.createPyPlot(axes, canvas.figure, newcplotB, 
                                    plotArgs['plot1'], **plotArgs)
    plotdataB = PyPlotB.data

    plotArgs['binname'] = 'luminosity_I'
    axes.setRanges(**plotArgs)
    axes.setProperties(**plotArgs)
    newcplotI = cplot.createCPlot(axes, **plotArgs['plot1'])  
    PyPlotI = pythplot.createPyPlot(axes, canvas.figure, newcplotI, 
                                    plotArgs['plot1'], **plotArgs)
    plotdataI = PyPlotI.data

    # convolving the luminosity grids with a gaussian to simulate seeing

    if gaussianFilter :
        gaussB = numpy.ndarray(shape=plotdataB.shape)
        gaussian_filter(plotdataB, 1, output=gaussB)#, mode="nearest")
        gaussI = numpy.ndarray(shape=plotdataI.shape)
        gaussian_filter(plotdataI, 1, output=gaussI)#, mode="nearest")
        for i in range(len(plotdataB)):
            for j in range(len(plotdataB[i])):
                plotdataB[i][j] = gaussB[i][j]
                pass
        for i in range(len(plotdataI)):
            for j in range(len(plotdataI[i])):
                plotdataI[i][j] = gaussI[i][j]
                pass

    # getting B-I color from luminosity grids

    for i in range(len(plotdataI)):
        for j in range(len(plotdataI[i])):
            try:
                plotdataI[i][j] = solLumB - 2.5*math.log10(plotdataB[i][j]) \
                    - (solLumI - 2.5*math.log10(plotdataI[i][j]))
                pass
            except:
                plotdataI[i][j] = 'nan'

    # plot the color data onto the axes
    PyPlotI.plotTheData()
    PyPlotI.colormap.set_bad('k')
    axes.addData(PyPlotI, plotArgs["plot1"])
    # get B-I color of total galaxy
    # M_B_tot = sumdata.getData('M_B')[snap-offset-1]
    # M_I_tot = sumdata.getData('M_I')[snap-offset-1]
    # BmI_tot = M_B_tot - M_I_tot

    # BmI_string = '$ (B-I)_{tot}={%5.2f}$'%BmI_tot
    # axes.text(0.05, 0.9, BmI_string, fontsize=13, transform=axes.transAxes, color='w')


    ### RENDER PLOT ###

    if renderPlot :

        plotArgs['cmin'], plotArgs['cmax'] = 0, 60
        plotArgs['plot1']['plotType']  = chyplot.CPyPlot.PLOT_RENDER
        plotArgs['binname'] = 'density'
        plotArgs['plot1']['nx']        = 30
        plotArgs['plot1']['ny']        = 30
        plotArgs['plot1']['type']      = 0
        plotArgs['plot1']['colormap']  = createColorBarFigs.getColorMapFromIndex(17)[1]
        plotArgs['plot1']['average']   = False
        plotArgs['plot1']['logscale']  = False
        plotArgs['plot1']['pointSize'] = 1
        plotArgs['plot1']['extent']    = 'smoothingLength'
        plotArgs['plot1']['normalize'] = True

        # plot rendered gas density
        PPlot.create(canvas.figure, renderAxes, **plotArgs )

        # plot stellar population
        for pop in range(0,nrOfPops):

            visitor = chyplot.cglobals.plmap.getSecond(constraints[pop])
            datacopy = data.limitsCopy(visitor, limits[pop][0], limits[pop][1], enums.T_star)
            # pdatacopy = glob.globals.dataManager.addData(datacopy)

            partData = [[0],[0]]

            partData[0] = datacopy.getDataArray(enums.T_star, 
                                                chyplot.cglobals.plmap.getSecond('x'), True)
            partData[1] = datacopy.getDataArray(enums.T_star, 
                                                chyplot.cglobals.plmap.getSecond('y'), True)

            try:
                renderAxes.scatter(partData[0], partData[1], c=colors[pop], s=size[pop])
            except:
                pass
            renderAxes.set_xlim(plotArgs["xmin"],plotArgs["xmax"])
            renderAxes.set_ylim(plotArgs["ymin"],plotArgs["ymax"])



    ### AESTHETICS ###

    # plot extra axes with time evolution / SFR
    if renderPlot :
        axes2 = canvas.figure.add_axes([0.22,0.05,0.585,0.12])
        axes2.plot(time,SF_data, color='0.4',zorder=0)
        # SFR_now = sum(sumdata.getData("SFR")[snap-offset-3:snap-offset+1])/4
        # time_now = sumdata.getData("time")[snap-offset-1]
        axes2.axvline(time_now,color='g')
        axes2.scatter(time_now, SFR_now, color=[0,1,0], s=10,zorder=10)
        axes2.text(0.92, 0.75, "$SFR$", fontsize=10.5, transform=axes2.transAxes, color='k')
        axes2.set_ylim(0,0.03)
        axes2.yaxis.set_major_locator(MaxNLocator(1))
        for ticklabel in axes2.get_yticklabels():
            ticklabel.set_visible(False)
    else :
        axes2 = canvas.figure.add_axes([0.145,0.13,0.585,0.12])
        axes2.plot(time,SF_data, color='0.4',zorder=0)
        # SFR_now = sum(sumdata.getData("SFR")[snap-offset-3:snap-offset+1])/4
        # time_now = sumdata.getData("time")[snap-offset-1]
        axes2.axvline(time_now,color='g')
        axes2.scatter(time_now, SFR_now, color=[0,1,0], s=10,zorder=15)
        axes2.text(0.85, 0.75, "$SFR$", fontsize=10.5, transform=axes2.transAxes, color='k')
        axes2.set_ylim(0,0.03)
        axes2.yaxis.set_major_locator(MaxNLocator(1))
        for ticklabel in axes2.get_xticklabels() + [axes2.get_yticklabels()[0]]:
            ticklabel.set_visible(False)
    
    # do various axes settings
    timestring = '$t={%5.2f}$'%data.time()
    if renderPlot :
        axes2.text(1.02, 0.1, timestring, fontsize=15, transform=axes2.transAxes, color='k')
    else :
        axes.text(0.8, 1.025, timestring, fontsize=fontsize, transform=axes.transAxes,color='k')
    axes.set_title("$B-I$", size=24)

    axes.set_aspect(1)

    axes.yaxis.set_major_locator(MaxNLocator(4))
    axes.xaxis.set_major_locator(MaxNLocator(4))
    
    canvas.figure.hideLabels(True)
    canvas.figure.subplots_adjust(left = 0.055, right = 0.97, bottom = 0.05, top = 0.98, wspace = 0, hspace = 0)

    # color bars
    if renderPlot :
        renderAxes.set_aspect(1)
        renderAxes.yaxis.set_major_locator(MaxNLocator(4))
        renderAxes.xaxis.set_major_locator(MaxNLocator(4))
        axes.hideColorbars()
        renderAxes.hideColorbars()
        
        # B-I plot
        cbplot = canvas.figure.getSubplot(1)
        cbdata, cbdict = cbplot.getData()[0]
        cbaxes = canvas.figure.add_axes([0.073,0.93,0.42,0.025])
        cbdata.updateColorBar(visible=True, cax=cbaxes, orientation='horizontal')
        #cbaxes.xaxis.set_major_locator(MaxNLocator(9))
        #cbaxes.xaxis.set_ticks([0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8])
        # render plot
        cbplot = canvas.figure.getSubplot(2)
        cbdata, cbdict = cbplot.getData()[0]
        cbaxes = canvas.figure.add_axes([0.533,0.93,0.42,0.025])
        cbdata.updateColorBar(visible=True, cax=cbaxes, orientation='horizontal')
        

    #cbaxes = canvas.figure.add_axes([0.91,0.05,0.02,0.91])
#cbaxes.yaxis.set_major_locator(LogLocator(base=10))
    #PyPlot.updateColorBar(visible=True, cax=cbaxes, orientation='vertical', ticks=[1,2,3,4,5], format='$%.f$')
#cbaxes.set_yscale('log')
#cbaxes.yaxis.set_minor_locator(LogLocator(base=10,subs=[2,3,4,5]))
    #cbaxes.set_ylabel('$stellar~ density~ \mathrm{[10^{6}M_{\odot}/kpc^{2}]}$', fontsize = 18)
#cbaxes.yaxis.set_major_locator(MaxNLocator(5))
#cbaxes.set_yticks([1])

    #axes.text(0.5,0.85,"$\mathrm{edge-on}$", fontsize=18, transform=axes.transAxes)
    #axes2.text(0.5,0.85,"$\mathrm{face-on}$", fontsize=18, transform=axes2.transAxes)

#canvas.figure.hideLabels(True)
#canvas.figure.subplots_adjust(left = 0.055, right = 0.91, bottom = 0.05, top = 0.96, wspace = 0, hspace = 0)

    name = str(run)+"_B-I_stellarbody_%04.f"%snap
    if gaussianFilter:
        name = name + "_gauss"+str(sigma)
    if renderPlot:
        name = name + "_render"
    name = name + ".png"
    if not os.path.exists("/home/robbert/Unief/Thesis/properties/figs/compare/" + str(run) + "/") :
    	    os.system('mkdir %s' %("/home/robbert/Unief/Thesis/properties/figs/compare/" + str(run) + "/"))
    figurepath = "/home/robbert/Unief/Thesis/properties/figs/compare/" + str(run) + "/" + name
    canvas.figure.savefig(figurepath, dpi=250, show=False)
    canvas.close()

#plt.show()







    # testing
#--------------------------------------------------------
#    for i in range(0,150):
#        for j in range(0,150):
#            if plotdataB[i][j] == 0.0:
#                #plotdataB[i][j] = 'nan'
#                pass
#            if plotdataI[i][j] == 0.0:
#                #plotdataI[i][j] = 'nan'
#                pass
    #plotdataB = numpy.zeros_like(plotdataB)
    #plotdataB[75][75] = 100
    #plotdataI = numpy.zeros_like(plotdataI)
    #axes.imshow(plotdataB)
    #plt.colorbar()
    #print plotdataB)
#    print plotdataB[0][0],plotdataB[75][75],plotdataI[0][0],plotdataI[75][75]
#    PyPlotB.data = gaussB
#    PyPlotI.data = gaussI
#--------------------------------------------------------

    # testing
    #print plotdataB[0][0],plotdataB[75][75],plotdataI[0][0],plotdataI[75][75]
    #axes.imshow(plotdataB)
    #PyPlotB.plotTheData()
    #PyPlotB.colormap.set_bad('k')
    #plt.show()
    #exit()

