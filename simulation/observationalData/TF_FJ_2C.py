import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from hyplot.visual.PRunData import PRunData, PFileData
from hyplot.plot import PFigure
import os
import chyplot
import enums
from pylab import *
from matplotlib.collections import LineCollection

matplotlib.rc('font', size=13)

sumdir = "/home/acloetos/Programs/scripts/sumsumfiles"
obsdir = "/home/acloetos/Programs/observational_data"
simdir = "/home/acloetos/Programs/gadget/output"
sim2dir = "/home/astro_dwarfgalaxies/acloetos/"
hyplotdir = "/home/acloetos/Programs/Hyplot/trunk/rundir/data/sumfiles"
inputsum = "/home/acloetos/Programs/scripts/inputsum"

isolated_10 = PFileData(os.path.join(sumdir, "sumsum_isolated_10"))
isolated_100 = PFileData(os.path.join(sumdir, "sumsum_isolated_100"))
#n50_e9 = PFileData(os.path.join(sumdir, "sumsum_Nsf50_Esn9_2C"))
m10_10_PI = PFileData(os.path.join(sumdir, "sumsum_mergers0.10E+10_PI"))
m25_10 = PFileData(os.path.join(sumdir, "sumsum_mergers0.25E+10_PI"))
m50_10 = PFileData(os.path.join(sumdir, "sumsum_mergers0.50E+10"))
m75_10 = PFileData(os.path.join(sumdir, "sumsum_mergers0.75E+10"))
m10_11_PI = PFileData(os.path.join(sumdir, "sumsum_mergers0.10E+11_PI"))
sims = [isolated_10, isolated_100, m10_10_PI,m25_10,m50_10,m75_10, m10_11_PI]
labels = ['$\mathrm{Isolated\ 10}$','$\mathrm{Isolated\ 100}$','$\mathrm{M_{ h}=1\ 10^{9}\ M_{\odot}}$','$\mathrm{M_{ h}=2.5\ 10^{9}\ M_{\odot}}$','$\mathrm{M_{ h}=5\ 10^{9}\ M_{\odot}}$','$\mathrm{M_{ h}=7.5\ 10^{9}\ M_{\odot}}$', '$\mathrm{M_{ h}=1\ 10^{10}\ M_{\odot}}$']
colors = ['r', 'b', 'Chartreuse','Chocolate', 'Orchid', 'Teal', 'Yellow']

input_isolated_10 = open(os.path.join(inputsum, "inputsum_isolated_10"))
input_isolated_100 = open(os.path.join(inputsum, "inputsum_isolated_100"))
#input_n50_e9 = open(os.path.join(inputsum, "inputsum_Nsf50_Esn9_2C"))
input_m10_10_PI = open(os.path.join(inputsum, "inputsum_0.10E+10_PI"))
input_m25_10 = open(os.path.join(inputsum, "inputsum_0.25E+10_PI"))
input_m50_10 = open(os.path.join(inputsum, "inputsum_0.50E+10"))
input_m75_10 = open(os.path.join(inputsum, "inputsum_0.75E+10"))
input_m10_11_PI = open(os.path.join(inputsum, "inputsum_0.10E+11_PI"))
input_sims = [input_isolated_10,input_isolated_100,input_m10_10_PI,input_m25_10,input_m50_10,input_m75_10, input_m10_11_PI]
lastsnap = []

for simul in input_sims:
    lastsnap.append([])
    for line in simul:
        temp=line.split('\t')
        if len(temp) > 2:
            if temp[1] != '0':
                lastsnap[-1].append(temp[1])

print "lastSnap: ", lastsnap            


TP00 = PFileData(os.path.join(obsdir, "tf_tullypierce.dat"))
K00 = PFileData(os.path.join(obsdir, "tf_kronawitter.dat"))
DR07 = PFileData(os.path.join(obsdir, "tf_vcsigma_dE.dat"))
MA01 = PFileData(os.path.join(obsdir, "tf_mag.dat"))
MG05 = PFileData(os.path.join(obsdir, "tf_mcgaugh.dat"))
VZ04tf = PFileData(os.path.join(obsdir, "tf_vanzee.dat"))
dI = PFileData(os.path.join(obsdir, "tf_dI.dat"))
DR05 = PFileData(os.path.join(obsdir, "derijcke2005.dat"))
GH03 = PFileData(os.path.join(obsdir, "geha2003.dat"))
GR03 = PFileData(os.path.join(obsdir, "graham2003.dat"))
KL05 = PFileData(os.path.join(obsdir, "kleyna2005.dat"))
MA98 = PFileData(os.path.join(obsdir, "mateo1998_modsigma.dat"))
PE93 = PFileData(os.path.join(obsdir, "Peterson1993.dat"))
WA07 = PFileData(os.path.join(obsdir, "walker2007_central_disp.dat"))
WI04 = PFileData(os.path.join(obsdir, "wilkinson2004_disp.dat"))
VZ04 = PFileData(os.path.join(obsdir, "vanZee2004.dat"))


G = 6.67300e-11 # in units of m^3 kg^-1 s^-2
kpc = 3.08568025e19 # conversion of kpc to m
solMass = 1.98892e36 # conversion of 10**6 solar mass to kg

visual.fig = plt.figure(FigureClass = PFigure.PFigure, figsize=(5,10))
ax1 = visual.fig.add_my_subplot(211)
ax2 = visual.fig.add_my_subplot(212)

ax1.set_ylabel('$\mathrm{\log_{10}(V_{c})\ [km/s]}$')
ax1.set_xlabel('$\mathrm{\log_{10}(L_{B})}$')
ax2.set_ylabel('$\mathrm{\log_{10}(\sigma)\ [km/s]}$')
ax2.set_xlabel('$\mathrm{\log_{10}(L_{B})}$')

ax1.text(0.9, 0.07, '$\mathrm{a.)}$', fontsize=14, transform=ax1.transAxes, color='k')
ax2.text(0.9, 0.07, '$\mathrm{b.)}$', fontsize=14, transform=ax2.transAxes, color='k')

ax1.plot(-0.4*(K00.getData('M_B')-5.48)-0.064, np.log10(K00.getData('Vcirc')),'.', markersize=5, color = 'y', label='$\mathrm{Observations}$', zorder=0)
ax1.plot(((-5.7611-15*dI.getData('dM')/670)-5.48)*(-0.4), -0.4715568862 + dI.getData('vc')*3/668, '.', markersize=5, color = 'y', zorder=0)
ax1.plot(MA01.getData('log10(LB)')+0.116, np.log10(MA01.getData('vc')), '.', markersize=5, color = 'y', zorder=0)
ax1.plot(np.log10(MG05.getData('Mstar')*1e10/MG05.getData('M/L')),np.log10(MG05.getData('Vcirc')), '.', markersize=5, color = 'y', zorder=0)
ax1.plot(-0.4*(TP00.getData('M_B')-5.48), TP00.getData('log10(Vcirc)')-0.30103, '.', markersize=5, color = 'y', zorder=0)
ax1.plot(-0.4*((VZ04tf.getData('d')-5*np.log10(16e5))-5.48), np.log10(np.sqrt((VZ04tf.getData('dv')*VZ04tf.getData('a')*0.1900698216)**2+VZ04tf.getData('sm')**2*2.45)) ,'.', markersize=5, color = 'y', zorder=0)
ax1.plot(np.log10(DR07.getData('Vcirc')), -0.4*(DR07.getData('M_B')-5.48), np.log10(DR07.getData('Vcirc')), '.', markersize=5, color = 'y', zorder=0)

print len(sims), len(colors), len(labels), len(lastsnap)

for sim, color, label, snaps in zip(sims, colors, labels, lastsnap):
    finalRun = []
    finalVc = []
    vcmaxList = []
    runs = sim.getData('run')
    LbList = sim.getData('log10(L_B)')
    print "sims: ", sim, color, label, snaps
    for run, snap in zip(runs, snaps):
        dr = chyplot.CDataGadget(run)
        if os.path.exists(os.path.join(simdir,"sim%04.d" %(run))):
            dr.setPrefix(os.path.join(simdir,"sim%04.d" %(run)))
        else:
            dr.setPrefix(os.path.join(sim2dir,"sim%04.d" %(run)))
        #dr.setPrefix(os.path.join(simdir,"sim%04.d" %(run)))
        #dr.checkFilesPresent()
        dr.set_file(int(snap))
        data2 = dr.readFile()
        chyplot.cglobals.plmap.setDataBlock(data2)
        data2.rcom(True, enums.T_dark)
        data2.vcom(True)
        data2.convertUnits()

        xmin = 0.1
        xmax = 21
        nbins = 1000
        binner = chyplot.CBinnerLinear()
        profile = chyplot.CProfile(xmin, xmax, nbins,
                                   "radius", "mass", enums.T_all)
        profile.computeProfile(data2, binner)
        profile.adaptiveBinning(50)
        bins_tot = profile.getBins()[:, 1] # y
        
        profile.cumulateProfile()
        cumul = profile.getProfile()
        vc = []
        for n in range(len(bins_tot)):
            vc.append(np.sqrt((G*cumul[n]*solMass)/(bins_tot[n]*kpc))/10**3)
        vcmaxList.append(max(vc))
        finalRun.append(int(run))
        finalVc.append(str(max(vc)))

    if len(vcmaxList) == len(LbList):
        if sim == sims[0] or sim == sims[1]:
            ax1.plot(LbList, np.log10(vcmaxList),'h--', color=color, markersize=5, label=label)
        else:
            ax1.plot(LbList, np.log10(vcmaxList),'h', color=color, markersize=5, label=label)

visual.fig.subplots_adjust(left = 0.10, right = 0.96, bottom = 0.11, top = 0.96,hspace = 0, wspace=0)

fontdict = {'size':11}
ax1.legend(loc=2, prop=fontdict, ncol=1, numpoints=1, scatterpoints=1, shadow=True, fancybox=True)

y2 = []
y3 = []
x = np.arange(1.5, 3., 0.1)
for i in range(len(x)):
    y2.append(3.42+3.09*x[i])
    y3.append(3.15+2.97*x[i])

ax1.plot(y2,x, '-', color='0.5')
ax1.plot(y3,x, '--', color='0.5')


ax1.set_ylim(1,2.8)
ax1.set_xlim(4.5, 11.5)

ax2.plot( ((GH03.getData('MV,0')+0.2752)/1.0345-5.441)/-2.5, np.log10(GH03.getData('disp')), '.', label='$\mathrm{Observations}$', markersize = 5, color='y', zorder=0)
ax2.plot((KL05.getData('M_B')-5.441)/-2.5, np.log10(KL05.getData('disp')), '.', markersize = 5, color='y', zorder=0)
ax2.plot((MA98.getData('M_B')-5.441)/-2.5, np.log10(MA98.getData('disp')), '.', markersize = 5, color = 'y', zorder=0)
ax2.plot((PE93.getData('M_B')-5.441)/-2.5, np.log10(PE93.getData('disp')), '.', markersize = 5, color = 'y', zorder=0)
ax2.plot( (-0.4*(VZ04.getData('mag_B') -5*np.log10(15.4e5)-5.48)), np.log10(VZ04.getData('sig')), '.', markersize = 5, color='y', zorder=0)
ax2.plot(DR05.getData('log10(L_B)'), np.log10(DR05.getData('<s>')), '.', markersize = 5, color='y', zorder=0)
#ax2.plot( n1_e1.getData('log10(L_B)'),n1_e1.getData('log10(sig)'),'h--',color='k', markersize=5, zorder=0, label='$\mathrm{n_{SF}=0.1,\ \epsilon_{FB}=0.1}$')
#ax2.plot( n6_e7.getData('log10(L_B)'),n6_e7.getData('log10(sig)'),'^--',color='b', markersize=5, zorder=0, label='$\mathrm{n_{SF}=6,\ \epsilon_{FB}=0.7}$')
#ax2.plot( n50_e9.getData('log10(L_B)'),n50_e9.getData('log10(sig)'),'d--',color='r', markersize=5, zorder=0, label='$\mathrm{n_{SF}=50,\ \epsilon_{FB}=0.9}$')

for sim, color, label in zip(sims, colors, labels):
    logLB = sim.getData('log10(L_B)')
    logSig = sim.getData('log10(sig)')
    x = []
    y = []
    if sim == sims[0] or sim == sims[1]:
        ax2.plot(logLB, logSig, 'h--', color=color, label=label)
    else:
        ax2.plot(logLB, logSig, 'h', color=color, label=label)

Lb = []
Lb_extended = []
x = np.arange(0, 3.1, 0.1)
y = np.arange(3, 12.1, 0.1)
for i in range(len(x)):
    Lb.append(6.02+1.57*x[i])
    Lb_extended.append(4.39+2.55*x[i])

#ax2.plot(x, Lb_extended)
ax2.set_ylim(0.7, 2.1)
ax2.set_xlim(4., 10.)


visual.fig.subplots_adjust(left = 0.13, right = 0.96, bottom = 0.06, top = 0.97,hspace = 0.15, wspace=0)
#visual.fig.hideLabels(True)

figurepath = "TF_FJ_2C.pdf"
visual.finalize(name=figurepath, dpi=200, show=False)
