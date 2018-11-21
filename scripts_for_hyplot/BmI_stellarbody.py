import os
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

matplotlib.rc('font', size=10)

solLumB = 5.441
solLumV = 4.820
solLumI = 4.148

edgeOn = False

### SETTINGS ###
run = '62002'
sumrun = run
snap_first = 85
snap_last  = 87
offset = 0
SFRmax = 0.025
maxSB = 28


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

# make plotargs dictionary

plotArgs = dict()

plotArgs['run']     = run
plotArgs['aspect']     = 'equal' 
if edgeOn:
        plotArgs['xname']      = 'x'
        plotArgs['yname']      = 'z'
        plotArgs['xtitle']     = 'x\ \mathrm{[kpc]}'
        plotArgs['ytitle']     = 'z\ \mathrm{[kpc]}'
else:
        plotArgs['xname']      = 'x'
        plotArgs['yname']      = 'y'
        plotArgs['xtitle']     = 'x\ \mathrm{[kpc]}'
        plotArgs['ytitle']     = 'y\ \mathrm{[kpc]}'

plotArgs['xmin'], plotArgs['xmax'] = -3,3
plotArgs['ymin'], plotArgs['ymax'] = -3,3
lengthSideBin = 0.02

plotArgs['plot1'] = dict()



gaussianFilter = True
sigma = 3

renderPlot = True

stellarPops = True
nrOfPops = 1      
constraints = ['age']#['age','age']#,'age']                             
limits = [[0,0.05]]#[[0.01,0.12],[0,0.01]]#,[0.05,0.1]]
colors = [(1,1,0)]#[(0.6,0.6,0),(1,1,0)]#,(0,0.6,0.6)]
size = [7]#[5,7]#,15]



# get analysis data from simulation

sumdata = PRunData(sumrun)
SF_data = sumdata.getData("SFR")
time = sumdata.getData("time")
M_B = sumdata.getData("M_B")
M_V = sumdata.getData("M_V")
M_I = sumdata.getData("M_I")



### PRODUCE PLOTS ###

print 'run',run
glob.globals.dataManager.setMaximumData( 10 )

for snap in range(snap_first, snap_last+1):


    ### DO SETUP ###

    # set up hyplot canvas
    if renderPlot:
        canvas = PCanvas.create(setupMDI=False, showCanvas=False, width = 3.5, height = 2.75)
    else :
        canvas = PCanvas.create(setupMDI=False, showCanvas=False, width = 5.5, height = 5.5)

    # read snapshot and set up data block
    dr = chyplot.CDataGadget(float(run))
    fdir = "/home/michele/sim/MoRIA/sim"+str(run)
    dr.setPrefix( fdir )
    dr.checkFilesPresent() # set the first and last dump
    dr.set_file( dr.firstDump() + snap)
    data = dr.readFile()
    print "time =", data.time()
    pdata = glob.globals.dataManager.addData(data)

    # do stuff to blocks
    data.rcom(True, enums.T_star, 0, 0, 0, True)
    data.vcom(True, enums.T_star)
    data.rotate(enums.T_gas, True, 10)
    data.convertUnits()
    chyplot.cglobals.plmap.setDataBlock(data)

    #data.rotate()
    
    # set up plotargs 
    plotArgs['fileNum'] = snap
    plotArgs['filename'] = glob.globals.dataManager.activePData().filename
    plotArgs['dataNumber'] = pdata.number

    # set up axes
    if renderPlot:
        axes   = canvas.figure.add_my_subplot(121)
        renderAxes = canvas.figure.add_my_subplot(122)
    else :
        axes   = canvas.figure.add_my_subplot(111)


    ### B-I COLOR PLOT ###

    plotArgs['cmin'], plotArgs['cmax'] = 0.6, 1.8 #1.5
    plotArgs['plot1']['plotType']  = chyplot.CPyPlot.PLOT_GRID
    plotArgs['plot1']['nx']        = int(math.ceil((plotArgs['xmax'] - plotArgs['xmin'])/lengthSideBin))
    plotArgs['plot1']['ny']        = int(math.ceil((plotArgs['ymax'] - plotArgs['ymin'])/lengthSideBin))
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
        gaussian_filter(plotdataB, sigma, output=gaussB)#, mode="nearest")
        gaussI = numpy.ndarray(shape=plotdataI.shape)
        gaussian_filter(plotdataI, sigma, output=gaussI)#, mode="nearest")
        for i in range(len(plotdataB)):
            for j in range(len(plotdataB[i])):
                plotdataB[i][j] = gaussB[i][j]
                pass
        for i in range(len(plotdataI)):
            for j in range(len(plotdataI[i])):
                plotdataI[i][j] = gaussI[i][j]
                pass

    # getting B-I color from luminosity grids
    
    for i in range(len(plotdataB)):
        for j in range(len(plotdataB[i])):
                if plotdataB[i][j] > 1e-9:
                        SB = solLumB + 21.572 - 2.5*math.log10(plotdataB[i][j]) + 15 + 2.5*math.log10(lengthSideBin**2)
                        #plotdataI[i][j] = SB
                        
                        SBI = solLumI + 21.572 - 2.5*math.log10(plotdataI[i][j]) + 15 + 2.5*math.log10(lengthSideBin**2)

                        if SBI < maxSB:
                                if SB < maxSB:
                                        plotdataI[i][j] = SB - SBI
                                else:
                                        plotdataI[i][j] = 'nan'#maxSB - SBI
                        else:
                                if SB < maxSB:
                                        plotdataI[i][j] =  SB - maxSB
                                else:
                                        plotdataI[i][j] = 'nan'
                else:
                        plotdataI[i][j] = 'nan'

    
                #else: 
                        #plotdataI[i][j] = 'nan'
    # plot the color data onto the axes
    PyPlotI.plotTheData()
    PyPlotI.colormap.set_bad('k')
    axes.addData(PyPlotI, plotArgs["plot1"])
    # get V-I color of total galaxy
    M_B_tot = M_B[snap-offset]
    #M_V_tot = M_V[snap-offset]
    M_I_tot = M_I[snap-offset]
    BmI_tot = M_B_tot - M_I_tot
    #VmI_tot = M_V_tot - M_I_tot

    BmI_string = '$ (B-I)_\mathrm{tot}={%5.2f}\ \mathrm{mag}$'%BmI_tot
    #VmI_string = '$ (V-I)_{tot}={%5.2f}$'%VmI_tot
    axes.text(0.05, 0.9, BmI_string, fontsize=10, transform=axes.transAxes, color='w')
    #axes.text(0.6, 0.9, VmI_string, fontsize=4, transform=axes.transAxes, color='w')


    ### RENDER PLOT ###

    if renderPlot :

        plotArgs['cmin'], plotArgs['cmax'] = 1e-26, 1e-22
        plotArgs['plot1']['plotType']  = chyplot.CPyPlot.PLOT_RENDER
        plotArgs['binname'] = 'density'
        plotArgs['plot1']['nx']        = 150
        plotArgs['plot1']['ny']        = 150
        plotArgs['plot1']['type']      = 0
        plotArgs['plot1']['colormap']  = createColorBarFigs.getColorMapFromIndex(21)[1]
        plotArgs['plot1']['average']   = False
        plotArgs['plot1']['logscale']  = True
        plotArgs['plot1']['pointSize'] = 1
        plotArgs['plot1']['extent']    = 'smoothingLength'
        plotArgs['plot1']['normalize'] = False

        # plot rendered gas density
        PPlot.create(canvas.figure, renderAxes, **plotArgs )

        visitorR = chyplot.cglobals.plmap.getSecond('radius')
        datacopyR =  data.limitsCopy(visitorR, 0, 5, enums.T_star)

        visitorT = chyplot.cglobals.plmap.getSecond('birthtime')
        datacopySFR =  datacopyR.limitsCopy(visitorT, data.time()-0.05, data.time(), enums.T_star)
        
        SFR_now= sum(datacopySFR.getDataArray(enums.T_star, chyplot.cglobals.plmap.getSecond('initialMass'), True))*1e6/(0.05*1e9)
        # plot stellar population
        
        for pop in range(0,nrOfPops):

            visitor = chyplot.cglobals.plmap.getSecond(constraints[pop])
            datacopy = data.limitsCopy(visitor, limits[pop][0], limits[pop][1], enums.T_star)
            pdatacopy = glob.globals.dataManager.addData(datacopy)

            partData = [[0],[0]]

            if edgeOn:
                    partData[0] = datacopy.getDataArray(enums.T_star, 
                                                        chyplot.cglobals.plmap.getSecond('x'), True)
                    partData[1] = datacopy.getDataArray(enums.T_star, 
                                                        chyplot.cglobals.plmap.getSecond('z'), True)

            else:
                    partData[0] = datacopy.getDataArray(enums.T_star, 
                                                        chyplot.cglobals.plmap.getSecond('x'), True)
                    partData[1] = datacopy.getDataArray(enums.T_star, 
                                                        chyplot.cglobals.plmap.getSecond('y'), True)

            try:
                renderAxes.scatter(partData[0], partData[1], c=colors[pop], s=size[pop], lw = 0.5)
            except:
                pass
            renderAxes.set_xlim(plotArgs["xmin"],plotArgs["xmax"])
            renderAxes.set_ylim(plotArgs["ymin"],plotArgs["ymax"])

        #u = numpy.arange(0.,2*math.pi, 0.01)
        #renderAxes.plot(numpy.cos(u), numpy.sin(u), '-')
        #axes.plot(numpy.cos(u), numpy.sin(u), '-')
        
    ### AESTHETICS ###

    # plot extra axes with time evolution / SFR

    SFR = sumdata.getData('SFR')
    time = sumdata.getData('time')
    SFR = [s*10**3 for s in SFR]
    averageOver = 1.
    SFRAverage = numpy.empty(int(math.ceil(len(SFR)/averageOver)))
    SFRAverage[0] = 0
    averageOver = int(averageOver)
    for i in xrange(int(math.ceil(len(SFR)/averageOver))-1):
        SFRAverage[i+1] = sum(SFR[averageOver*(i+1)-averageOver+1:averageOver*(i+1)+1])/averageOver
    timeAverage = time[::averageOver]
    
    fontsize = 10

    if renderPlot :
        axes2 = canvas.figure.add_axes([0.127,0.05,0.825,0.15])
        axes2.plot(timeAverage,SFRAverage, color='k',zorder=0, linewidth=1)
        #SFR_now = sumdata.getData("SFR")[snap-offset]
        time_now = sumdata.getData("time")[snap-offset]
        axes2.axvline(time_now,color='g')
        axes2.tick_params(axis='both', which='major', labelsize=5)
        #axes2.scatter(time_now, SFR_now, color=[0,1,0], s=10,zorder=10)
        axes2.text(0.05, 0.5, "$SFR =$ " + str(round(SFR_now,3)) + " $M_\odot$/yr", fontsize=fontsize, transform=axes2.transAxes, color='k')
        axes2.set_ylim(0,SFRmax*10**3)
        axes2.set_xlim(0,sumdata.getData('time')[-1])
        axes2.set_xlabel('$t\ \mathrm{[Gyr]}$', fontsize=fontsize)
        axes2.set_ylabel('$\mathrm{SFR}\ [10^{-3}\mathrm{M}_\odot\mathrm{/yr]}$', fontsize=fontsize)
        axes2.yaxis.set_major_locator(MaxNLocator(1))
        #for ticklabel in axes2.get_yticklabels():
            #ticklabel.set_visible(False)
    else :
        axes2 = canvas.figure.add_axes([0.15,0.05,0.8,0.08])
        axes2.plot(time,SF_data, color='0.4',zorder=0)
        SFR_now = sumdata.getData("SFR")[snap-offset]
        time_now = sumdata.getData("time")[snap-offset]

        axes2.axvline(time_now,color='g')
        axes2.scatter(time_now, SFR_now, color=[0,1,0], s=10,zorder=10)
        axes2.text(0.05, 0.75, "$SFR =$ " + str(round(SFR_now,3)) + " $M_\odot$/yr", fontsize=fontsize, transform=axes2.transAxes, color='k')
        axes2.set_ylim(0,SFRmax)
        axes2.set_xlim(sumdata.getData('time')[0],sumdata.getData('time')[-1])
        axes2.yaxis.set_major_locator(MaxNLocator(1))
        for ticklabel in axes2.get_yticklabels():
            ticklabel.set_visible(False)
    
    # do various axes settings
    timestring = '$t={%5.2f}\ \mathrm{Gyr}$'%data.time()
    if renderPlot :
        axes2.text(0.37, 1.2, timestring, fontsize=7., transform=axes2.transAxes, color='k')
    else :
        props = dict(boxstyle='round', facecolor='w', alpha=0.9)
        axes.text(0.7, 0.9, timestring, fontsize=7, transform=axes.transAxes, color='k', bbox=props)
    axes.set_title("$B-I$", size=24)

    axes.set_aspect(1)

    axes.yaxis.set_major_locator(MaxNLocator(4))
    axes.xaxis.set_major_locator(MaxNLocator(4))
    
    canvas.figure.hideLabels(True)
    canvas.figure.subplots_adjust(left = 0.13, right = 0.95, bottom = 0.21, top = 0.97, wspace = 0, hspace = 0)

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
        cbaxes = canvas.figure.add_axes([0.127,0.905,0.408,0.025])
        cbdata.updateColorBar(visible=True, cax=cbaxes, orientation='horizontal')
        cbaxes.tick_params(axis='both', which='major', labelsize=4)
        cbaxes.text(0.25,1.85, '$B-I\ \mathrm{[mag/arcsec}^2]$', fontsize=6, transform=cbaxes.transAxes, color='k')
        #cbaxes.xaxis.set_ticks([0.6,0.9,1.2,1.5,1.8])
        #cbaxes.xaxis.set_major_locator(MaxNLocator(9))
        #cbaxes.xaxis.set_ticks([0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8])
        # render plot
        cbplot = canvas.figure.getSubplot(2)
        cbdata, cbdict = cbplot.getData()[0]
        cbaxes = canvas.figure.add_axes([0.542,0.905,0.408,0.025])
        cbaxes.tick_params(axis='both', which='major', labelsize=4)


        cbdata.updateColorBar(visible=True, cax=cbaxes, orientation='horizontal')
        cbaxes.text(0.2,1.55, '$\Sigma_\mathrm{gas}\ [10^6\mathrm{M}_\odot\mathrm{/kpc}^2]$', fontsize=6, transform=cbaxes.transAxes, color='k')

    else:
        axes.hideColorbars()
        cbplot = canvas.figure.getSubplot(1)
        cbdata, cbdict = cbplot.getData()[0]
        cbaxes = canvas.figure.add_axes([0.02,0.3,0.025,0.42])
        cbdata.updateColorBar(visible=True, cax=cbaxes, orientation='vertical')
        cbaxes.xaxis.set_major_locator(MaxNLocator(9))
        cbaxes.xaxis.set_ticks([0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8])


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
    if renderPlot:
        name += "_render"
    name = name + ".png"
    if edgeOn:
        directory = "/home/michele/sim/analysis/results/colorMapsXZ/sim" + str(run) + "/"
    else:
        directory = "/home/michele/sim/analysis/results/colorMaps/sim" + str(run) + "/"
    if not os.path.exists(directory) :
        os.system('mkdir -p %s' %(directory))
    figurepath = directory + name
    canvas.figure.savefig(figurepath, dpi=600, show=False)
    # canvas.close()

    plt.show()







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

