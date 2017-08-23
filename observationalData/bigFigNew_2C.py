import numpy as np
from hyplot.visual.PRunData import PRunData, PFileData
from hyplot.plot import PFigure
import matplotlib
import matplotlib.pyplot as plt
from pylab import *
from matplotlib.collections import LineCollection

matplotlib.rc('font', size=14)

sumdir = "/home/acloetos/Programs/scripts/sumsumfiles"
obsdir = "/home/acloetos/Programs/observational_data"
sumsumdir = "/home/acloetos/Programs/Hyplot/trunk/rundir/data/sumfiles"

n1_e1 = PFileData(os.path.join(sumdir, "sumsum_Nsf1_Esn1_2C"))
n6_e1 = PFileData(os.path.join(sumdir, "sumsum_Nsf6_Esn1_2C"))
n6_e3 = PFileData(os.path.join(sumdir, "sumsum_Nsf6_Esn3_2C"))
n6_e5 = PFileData(os.path.join(sumdir, "sumsum_Nsf6_Esn5_2C"))
n6_e7 = PFileData(os.path.join(sumdir, "sumsum_Nsf6_Esn7_2C"))
n6_e9 = PFileData(os.path.join(sumdir, "sumsum_Nsf6_Esn9_2C"))
n50_e3 = PFileData(os.path.join(sumdir, "sumsum_Nsf50_Esn3_2C"))
n50_e5 = PFileData(os.path.join(sumdir, "sumsum_Nsf50_Esn5_2C"))
n50_e7 = PFileData(os.path.join(sumdir, "sumsum_Nsf50_Esn7_2C"))
n50_e9 = PFileData(os.path.join(sumdir, "sumsum_Nsf50_Esn9_2C"))

DR05 = PFileData(os.path.join(obsdir , "vsigeps.dat"))
GG03 = PFileData(os.path.join(obsdir , "graham_LB_Re.dat"))
GG03extraextra = PFileData(os.path.join(obsdir, "graham2003.dat"))
GG03extra = PFileData(os.path.join(obsdir, "Grahamfig8mumb.dat"))
MW = PFileData(os.path.join(obsdir, "lgdsph_MW.dat"))
M31 = PFileData(os.path.join(obsdir, "lgdsph_M31.dat"))
Free = PFileData(os.path.join(obsdir, "lgdsph_free.dat"))
Perseus = PFileData(os.path.join(obsdir, "perseus2.dat")) 
A = PFileData(os.path.join(obsdir, "antlia.dat")) 
MA98 = PFileData(os.path.join(obsdir , "mateo1998_modsigma.dat"))
GE03 = PFileData(os.path.join(obsdir , "geha2003.dat"))
PE03 = PFileData(os.path.join(obsdir , "Peterson1993.dat"))
MA98 = PFileData(os.path.join(obsdir , "mateo1998_modsigma.dat"))
NA04 = PFileData(os.path.join(obsdir, "naga2004.dat"))
MI07 = PFileData(os.path.join(obsdir, "michielsen2007.dat"))
gr   = PFileData(os.path.join(obsdir, "grebel2003_clean.dat"))
dsph = PFileData(os.path.join(obsdir, "grebel2003_dSph.dat"))
de   = PFileData(os.path.join(obsdir, "grebel2003_dE.dat"))
dirr = PFileData(os.path.join(obsdir, "grebel2003_dIrr.dat"))
dirrdsph = PFileData(os.path.join(obsdir, "grebel2003_dirrdsph.dat"))
dsphShar = PFileData(os.path.join(obsdir, "sharina2008_dSph.dat"))
dirrShar = PFileData(os.path.join(obsdir, "sharina2008_dIrr.dat"))
dirrdsphShar = PFileData(os.path.join(obsdir, "sharina2008_dSphdIrr.dat"))
dsphLian = PFileData(os.path.join(obsdir, "lianou2010_dSph.dat"))
dirrdsphLian = PFileData(os.path.join(obsdir, "lianou2010_dSphdIrr.dat"))

dataBU97 = PFileData(os.path.join(obsdir, "burstein1997.dat.new"))
dataDR05 = PFileData(os.path.join(obsdir, "derijcke2005.dat"))

visual.fig = plt.figure(FigureClass = PFigure.PFigure, figsize=(10,12))
ax1 = visual.fig.add_my_subplot(321)
ax2 = visual.fig.add_my_subplot(322)
ax3 = visual.fig.add_my_subplot(323, sharex=ax1)
ax4 = visual.fig.add_my_subplot(324)
ax5 = visual.fig.add_my_subplot(325, sharex=ax1)
ax6 = visual.fig.add_my_subplot(326, sharex=ax4)

ax1.text(0.9, 0.07, '$\mathrm{a.)}$', fontsize=14, transform=ax1.transAxes, color='k')
ax2.text(0.9, 0.07, '$\mathrm{b.)}$', fontsize=14, transform=ax2.transAxes, color='k')
ax3.text(0.9, 0.07, '$\mathrm{c.)}$', fontsize=14, transform=ax3.transAxes, color='k')
ax4.text(0.9, 0.07, '$\mathrm{d.)}$', fontsize=14, transform=ax4.transAxes, color='k')
ax5.text(0.9, 0.07, '$\mathrm{e.)}$', fontsize=14, transform=ax5.transAxes, color='k')
ax6.text(0.9, 0.07, '$\mathrm{f.)}$', fontsize=14, transform=ax6.transAxes, color='k')

## Mv_Re
#observational data 
ax1.plot(DR05.getData('MB')-0.7, np.log10(DR05.getData('re_in_kpc')),  '.y', markersize=6, label='$\mathrm{Observations}$' , zorder=0)
ax1.plot(1.0345*GG03.getData('M_B')-0.2752, GG03.getData('log10(Re)'), '.y', markersize=6, zorder=0)
ax1.plot(MW.getData('M_V'), np.log10(np.tan(MW.getData("re\"")/60.*np.pi/180.)*MW.getData("D")), '.', c='y', markersize=6, zorder=0)
ax1.plot(M31.getData('M_V'), np.log10(M31.getData('r_e_L')), '.y', markersize=6, zorder=0)
ax1.plot(Free.getData('M_V'), np.log10(Free.getData('r_e_L')), '.y', markersize=6, zorder=0)
ax1.plot(Perseus.getData('m555') - 0.036 - 0.051*(Perseus.getData('m555_1re') - Perseus.getData('m814_1re')) - 34.26, np.log10(Perseus.getData('Re')), '.y', markersize=6, zorder=0)
ax1.plot(A.getData('dm')-2.63*A.getData('ebv')+0.183*A.getData('dct')+0.208-32.73, np.log10(350*np.tan(A.getData('dre')/3600.*np.pi/180.)),'.y', markersize=6, zorder=0) 
ax1.plot(n1_e1.getData('M_B')-n1_e1.getData('M_B-M_V'), n1_e1.getData('log10(r_e_L)'),  'h--', c='k', markersize=5, label='$\mathrm{n_{SF}=0.1\ cm^{-3},\ \epsilon_{SF}=0.1}$' , zorder=0)

#fontdict = {'size':12}
#ax1.legend(loc=2, prop=fontdict, ncol=1, numpoints=1, scatterpoints = 1, shadow=True, fancybox=True)

#final
MV6 = [n6_e1.getData('M_V'),n6_e3.getData('M_V'),n6_e5.getData('M_V'),n6_e7.getData('M_V'),n6_e9.getData('M_V')]
Mb6 = [n6_e1.getData('M_B'),n6_e3.getData('M_B'),n6_e5.getData('M_B'),n6_e7.getData('M_B'),n6_e9.getData('M_B')]
Mb50 = [n50_e3.getData('M_B'),n50_e5.getData('M_B'),n50_e7.getData('M_B'),n50_e9.getData('M_B')]
LogRe6 = [n6_e1.getData('log10(r_e_L)'), n6_e3.getData('log10(r_e_L)'), n6_e5.getData('log10(r_e_L)'), n6_e7.getData('log10(r_e_L)'), n6_e9.getData('log10(r_e_L)')]
MV50 = [n50_e3.getData('M_V'),n50_e5.getData('M_V'),n50_e7.getData('M_V'),n50_e9.getData('M_V')]
LogRe50 = [n50_e3.getData('log10(r_e_L)'), n50_e5.getData('log10(r_e_L)'), n50_e7.getData('log10(r_e_L)'), n50_e9.getData('log10(r_e_L)')]
Re6 = [n6_e1.getData('r_e_L'),n6_e3.getData('r_e_L'),n6_e5.getData('r_e_L'),n6_e7.getData('r_e_L'),n6_e9.getData('r_e_L')]
Re50 = [n50_e3.getData('r_e_L'),n50_e5.getData('r_e_L'),n50_e7.getData('r_e_L'),n50_e9.getData('r_e_L')]
Ie6 = [n6_e1.getData('Ie'),n6_e3.getData('Ie'),n6_e5.getData('Ie'),n6_e7.getData('Ie'),n6_e9.getData('Ie')]
Ie50 = [n50_e3.getData('Ie'),n50_e5.getData('Ie'),n50_e7.getData('Ie'),n50_e9.getData('Ie')]
Sigma6 = [n6_e1.getData('sig_1D_c'),n6_e3.getData('sig_1D_c'),n6_e5.getData('sig_1D_c'),n6_e7.getData('sig_1D_c'),n6_e9.getData('sig_1D_c')]
Sigma50 = [n50_e3.getData('sig_1D_c'),n50_e5.getData('sig_1D_c'),n50_e7.getData('sig_1D_c'),n50_e9.getData('sig_1D_c')]
Lb6 = [n6_e1.getData('L_B'),n6_e3.getData('L_B'),n6_e5.getData('L_B'),n6_e7.getData('L_B'),n6_e9.getData('L_B')]
Lb50 = [n50_e3.getData('L_B'),n50_e5.getData('L_B'),n50_e7.getData('L_B'),n50_e9.getData('L_B')]
VminI6 = [n6_e1.getData('V-I'),n6_e3.getData('V-I'),n6_e5.getData('V-I'),n6_e7.getData('V-I'),n6_e9.getData('V-I')]
VminI50 = [n50_e3.getData('V-I'),n50_e5.getData('V-I'),n50_e7.getData('V-I'),n50_e9.getData('V-I')]
mu06 = [n6_e1.getData('mu0'),n6_e3.getData('mu0'),n6_e5.getData('mu0'),n6_e7.getData('mu0'),n6_e9.getData('mu0')]
mu050 = [n50_e3.getData('mu0'),n50_e5.getData('mu0'),n50_e7.getData('mu0'),n50_e9.getData('mu0')]
n6 = [n6_e1.getData('n'),n6_e3.getData('n'),n6_e5.getData('n'),n6_e7.getData('n'),n6_e9.getData('n')]
n50 = [n50_e3.getData('n'),n50_e5.getData('n'),n50_e7.getData('n'),n50_e9.getData('n')]
FeH6 =[n6_e1.getData('[Fe/H]'), n6_e3.getData('[Fe/H]'), n6_e5.getData('[Fe/H]'), n6_e7.getData('[Fe/H]'), n6_e9.getData('[Fe/H]')]
FeH50 = [n50_e3.getData('[Fe/H]'), n50_e5.getData('[Fe/H]'), n50_e7.getData('[Fe/H]'), n50_e9.getData('[Fe/H]')]

conc1 = []
x = []
y = []
c = []
telmod1 = xrange(5)
tempList1 = [0.1,0.3,0.5,0.7,0.9]
for l in telmod1:
    Mv6List = []
    LogRe6List = []
    for i in range(len(MV6[l])):
        Mv6List.append([MV6[l][i]])
        LogRe6List.append([LogRe6[l][i]])
        x.append(MV6[l][i])
        y.append(LogRe6[l][i])
        c.append(tempList1[l])
    conc1.append(np.append(Mv6List, LogRe6List, axis=1))
linesN6 = LineCollection([conc1[l] for l in telmod1],linewidth=1,cmap=cm.winter)
res6 = ax1.scatter(x,y,s=20,c=c,facecolor= 'w', marker='d',cmap=cm.winter, linewidths=0.5, label = '$\mathrm{n_{SF}=6\ cm^{-3}}$')

telmod2 = xrange(4)
conc2 = []
x = []
y = []
c = []

tempList2 = [0.3,0.5,0.7,0.9]
for k in telmod2:
    Mv50List = []
    LogRe50List = []
    for m in range(len(Mb50[k])):
        Mv50List.append([MV50[k][m]])
        LogRe50List.append([LogRe50[k][m]])
        x.append(MV50[k][m])
        y.append(LogRe50[k][m])
        c.append(tempList2[k])
    conc2.append(np.append(Mv50List, LogRe50List, axis=1))
linesN50 = LineCollection([conc2[k] for k in telmod2],linewidth=1,cmap=cm.autumn)
res50 = ax1.scatter(x,y,s=20,c=c,marker='^', facecolor='w', cmap=cm.autumn, linewidths=0.5, label='$\mathrm{n_{SF}=50\ cm^{-3}}$')

fontdict = {'size':12}
ax1.legend(loc=2, prop=fontdict, ncol=1, numpoints=1, scatterpoints = 1, shadow=True, fancybox=True)

epsilon1=[0.1,0.3,0.5,0.7,0.9]
arr = np.array([epsilon1[l] for l in telmod1])
linesN6.set_array(arr)
ax1.add_collection(linesN6)

ax6.set_position([.53,.60,0.45,0.25])

cb6 = visual.fig.colorbar(res6, orientation = 'vertical', shrink=1, ticks = [0.1,0.3,0.5,0.7,0.9])
cb6.set_label('$\mathrm{\epsilon_{\\rm FB}~ (n_{\\rm SF}=6~ cm^{-3})}$',size=14)

epsilon2=[0.3,0.5,0.7,0.9]
arr = np.array([epsilon2[k] for k in telmod2])
linesN50.set_array(arr)
ax1.add_collection(linesN50)

ax6.set_position([.53,.18,0.45,0.25])
cb50 = visual.fig.colorbar(res50, orientation='vertical', ticks=[0.3,0.5,0.7,0.9], pad=0.10, shrink=1)
cb50.set_label('$\mathrm{\epsilon_{\\rm FB}~ (n_{\\rm SF}=50~ cm^{-3})}$',size=14)

##resFP
#definitions
def FP_nolog(Ie, sigma_c):
    return -0.629 - 0.845 * Ie + 1.379*sigma_c

def FP(Ie, sigma_c):
    global np, FP_nolog
    return FP_nolog(np.log10(Ie), np.log10(sigma_c))

def FPdev(Ie, sigma_c, re):
    global np, FP_nolog
    return FP_nolog(np.log10(Ie), np.log10(sigma_c)) - np.log10(re)

def FPdev_nolog(Ie, sigma_c, re):
    return -0.629 - 0.845 * Ie + 1.379*sigma_c - re
    
ax2.set_xlabel('$\mathrm{\log_{10}(L_{B})}$')
ax2.set_ylabel('$\mathrm{FP - \log_{10}(R_{e})}$')

ax2.plot((dataBU97.getData('15')[:400]), FPdev_nolog(dataBU97.getData('log10(Ie)')[:400],dataBU97.getData('log10(sigma_c)')[:400], dataBU97.getData('log10(Re)')[:400]), '.', c='y', markersize=6,label='$\mathrm{Observations}$')
ax2.plot((dataDR05.getData('log10(L_B)')), FPdev(dataDR05.getData('Ie'),dataDR05.getData('<s>'), dataDR05.getData('r_e_L')), '.', c='y', markersize=6 )

ax2.plot(np.log10(n1_e1.getData('L_B')), FPdev(n1_e1.getData('Ie'),n1_e1.getData('sig_1D_c'), n1_e1.getData('r_e_L')), 'h--',color='k',  markersize=4,label='$\mathrm{n_{SF}=0.1\,\ \epsilon_{FB}=0.1}$' )

## FPdev N=6 data ##
conc1 = []
x = []
y = []
c = []
for l in telmod1:
    Lb6List = []
    FPdev6List = []
    for i in range(len(Lb6[l])):
        x.append(np.log10(Lb6[l][i]))
        y.append(FPdev(Ie6[l][i],Sigma6[l][i],Re6[l][i]))
        Lb6List.append([np.log10(Lb6[l][i])])
        FPdev6List.append([FPdev(Ie6[l][i],Sigma6[l][i],Re6[l][i])])
        c.append(tempList1[l])
    conc1.append(np.append(Lb6List, FPdev6List, axis=1))
linesN6 = LineCollection([conc1[l] for l in telmod1],linewidth=1,cmap=cm.winter)
res6 = ax2.scatter(x,y,s=20,c=c,marker='d',cmap=cm.winter,facecolor= 'w', linewidths=0.5, label = '$\mathrm{n_{SF}=6\ cm^{-3}}$')

## FPdev N=50 data ##
conc2 = []
x = []
y = []
c = []
for k in telmod2:
    Lb50List = []
    FPdev50List = []
    for m in range(len(Lb50[k])):
        x.append(np.log10(Lb50[k][m]))
        y.append(FPdev(Ie50[k][m],Sigma50[k][m],Re50[k][m]))
        Lb50List.append([np.log10(Lb50[k][m])])
        FPdev50List.append([FPdev(Ie50[k][m],Sigma50[k][m],Re50[k][m])])
        c.append(tempList2[k])
    conc2.append(np.append(Lb50List, FPdev50List, axis=1))
linesN50 = LineCollection([conc2[k] for k in telmod2],linewidth=1,cmap=cm.autumn)
res50 = ax2.scatter(x,y,s=20,c=c,marker='^',cmap=cm.autumn,facecolor= 'w', linewidths=0.5, label = '$\mathrm{n_{SF}=50\ cm^{-3}}$')

x = np.arange(3, 14, 1)
y = np.zeros(len(x))
ax2.plot(x, y, 'k-')

epsilon1=[0.1,0.3,0.5,0.7,0.9]
arr = np.array([epsilon1[l] for l in telmod1])
linesN6.set_array(arr)
ax2.add_collection(linesN6)

epsilon2=[0.3,0.5,0.7,0.9]
arr = np.array([epsilon2[k] for k in telmod2])
ax2.add_collection(linesN50)
linesN50.set_array(arr)

ax2.set_xlim(4, 12)
ax2.set_ylim(-0.7, 1.5)

## MV_VminI
#observational data
ax3.plot(DR05.getData('MB')-0.7, 0.198*DR05.getData('dfe') + 1.207, '.', c='y', markersize=6, zorder=0)
ax3.plot(MW.getData('MV'), 0.198*MW.getData('[Fe/H]') + 1.207, '.', c='y', markersize=6, zorder=0)
ax3.plot(M31.getData('MV'),M31.getData('mu_0'), '.', c='y', markersize=6, zorder=0)
ax3.plot(Free.getData('MV'),Free.getData('mu_0'), '.', c='y', markersize=6, zorder=0)
ax3.plot(Perseus.getData('m555') - 0.036 - 0.051*(Perseus.getData('m555_1re') - Perseus.getData('m814_1re')) - 34.26, Perseus.getData('m555_1re') - 0.036 - 0.051*(Perseus.getData('m555_1re') - Perseus.getData('m814_1re')) - 34.26 - (Perseus.getData('m814_1re') - 0.443 - 0.002*(Perseus.getData('m555_1re') - Perseus.getData('m814_1re')) - 34.26), '.', c='y', markersize=6, zorder=0)
ax3.plot(A.getData('dm')-2.63*A.getData('ebv')+0.183*A.getData('dct')+0.208-32.73, 0.465*(A.getData('dct')-1.97*A.getData('ebv')) + 0.338, '.', c='y', markersize=6, zorder=0)
ax3.plot(n1_e1.getData('M_V'), n1_e1.getData('V-I'),  'h--', c='k', markersize=5, label='$\mathrm{n_{SF}=0.1\ cm^{-3},\ \epsilon_{SF}=0.1}$', zorder=0)

conc1 = []
x = []
y = []
c = []
for l in telmod1:
    Mv6List = []
    VminI6List = []
    for i in range(len(MV6[l])):
        Mv6List.append([MV6[l][i]])
        VminI6List.append([VminI6[l][i]])
        x.append(MV6[l][i])
        y.append(VminI6[l][i])
        c.append(tempList1[l])
    conc1.append(np.append(Mv6List, VminI6List, axis=1))
linesN6 = LineCollection([conc1[l] for l in telmod1],linewidth=1,cmap=cm.winter)
res6 = ax3.scatter(x,y,s=20,c=c,facecolor= 'w', marker='d',cmap=cm.winter, linewidths=0.5, label = '$\mathrm{n_{SF}=6\ cm^{-3}}$')

epsilon1=[0.1,0.3,0.5,0.7,0.9]
arr = np.array([epsilon1[l] for l in telmod1])
linesN6.set_array(arr)
ax3.add_collection(linesN6)

conc2 = []
x = []
y = []
c = []
for k in telmod2:
    Mv50List = []
    VminI50List = []
    for m in range(len(Mb50[k])):
        Mv50List.append([MV50[k][m]])
        VminI50List.append([VminI50[k][m]])
        x.append(MV50[k][m])
        y.append(VminI50[k][m])
        c.append(tempList2[k])
    conc2.append(np.append(Mv50List, VminI50List, axis=1))
linesN50 = LineCollection([conc2[k] for k in telmod2],linewidth=1,cmap=cm.autumn)
res50 = ax3.scatter(x,y,s=20,c=c,marker='^', facecolor='w', cmap=cm.autumn, linewidths=0.5, label='$\mathrm{n_{SF}=50\ cm^{-3}}$')

epsilon2=[0.3,0.5,0.7,0.9]
arr = np.array([epsilon2[k] for k in telmod2])
linesN50.set_array(arr)
ax3.add_collection(linesN50)

## MV_mu0
#observatinal data
def f1(MB):
    return 1.0345*MB-0.2752
def f2(mu0, MB):
    return mu0+0.0345*MB-0.2752

ax6.plot(DR05.getData('M_B')-0.7, DR05.getData('m0')+0.55,  '.y', markersize=6, label='$\mathrm{DR05}$', zorder=0)
ax6.plot(f1(GG03extra.getData('M_B')), f2(GG03extra.getData('mu_0'),GG03extra.getData('M_B')), '.y', markersize=6, zorder=0)
ax6.plot(MW.getData('MV'), MW.getData('sb0'),'.y', markersize=6, zorder=0)
ax6.plot(M31.getData('MV'), M31.getData('sb0'),'.y', markersize=6, zorder=0)
ax6.plot(Free.getData('MV'), Free.getData('sb0'),'.y', markersize=6, zorder=0)
ax6.plot(Perseus.getData('m555') - 0.036 - 0.051*(Perseus.getData('m555_1re') - Perseus.getData('m814_1re')) - 34.26, Perseus.getData('dmu') - 0.036 - 0.051*(Perseus.getData('m555_1re')-Perseus.getData('m814_1re')) -0.54,".y", markersize=6, zorder=0)

x = []
y = []
c = []
conc1 = []
for l in telmod1:
    Mv6List = []
    Mu06List = []
    for i in range(len(MV6[l])):
        Mv6List.append([MV6[l][i]])
        Mu06List.append([mu06[l][i]])
        x.append(MV6[l][i])
        y.append(mu06[l][i])
        c.append(tempList1[l])
    conc1.append(np.append(Mv6List, Mu06List, axis=1))
linesN6 = LineCollection([conc1[l] for l in telmod1],linewidth=1,cmap=cm.winter)
res6 = ax6.scatter(x,y,s=20,c=c,marker='d',cmap=cm.winter, linewidths=0.5)

epsilon1=[0.1,0.3,0.5,0.7,0.9]
arr = np.array([epsilon1[l] for l in telmod1])
linesN6.set_array(arr)
ax6.add_collection(linesN6)

x = []
y = []
c = []
conc2 = []
for k in telmod2:
    Mv50List = []
    Mu050List = []
    for m in range(len(MV50[k])):
        Mv50List.append([MV50[k][m]])
        Mu050List.append([mu050[k][m]])
        x.append(MV50[k][m])
        y.append(mu050[k][m])
        c.append(tempList2[k])
    conc2.append(np.append(Mv50List, Mu050List, axis=1))
linesN50 = LineCollection([conc2[k] for k in telmod2],linewidth=1,cmap=cm.autumn)
res50 = ax6.scatter(x,y,s=20,c=c,marker='^',cmap=cm.autumn, linewidths=0.5)

epsilon2=[0.3,0.5,0.7,0.9]
arr = np.array([epsilon2[k] for k in telmod2])
linesN50.set_array(arr)
ax6.add_collection(linesN50)

ax6.plot(-n1_e1.getData('M_B-M_V')+n1_e1.getData('M_B'), n1_e1.getData('mu0'), 'h--', c='k', markersize=5)

## Mv_n
# Observational data
ax5.plot(DR05.getData('M_B')-0.7, DR05.getData('n')+0.55,  '.', c='y', markersize=6, label='$\mathrm{Observations}$', zorder=0)
ax5.plot(1.0345*GG03extraextra.getData('M_B')-0.2752 , GG03extraextra.getData('n'), '.', c='y', markersize=6, zorder=0)
ax5.plot(MW.getData('MV'), MW.getData('n'),'.', c='y', markersize=6, zorder=0)
ax5.plot(M31.getData('MV'), M31.getData('n'),'.', c='y', markersize=6, zorder=0)
ax5.plot(Free.getData('MV'), Free.getData('n'),'.', c='y', markersize=6, zorder=0)
ax5.plot(Perseus.getData('m555') - 0.036 - 0.051*(Perseus.getData('m555_1re') - Perseus.getData('m814_1re')) - 34.26, Perseus.getData('n'), ".y", markersize=6, zorder=0)

ax5.plot(-n1_e1.getData('M_B-M_V')+n1_e1.getData('M_B'), n1_e1.getData('n'), 'h--', c='k', markersize=4, label='$\mathrm{n_{SF}= 0.1,\ \epsilon_{FB}=0.1}$')

# Simulation data
def f1(MB):
    return 1.0345*MB-0.2752
def f2(mu0, MB):
    return mu0+0.0345*MB-0.2752

conc1 = []
x = []
y = []
c = []
for l in telmod1:
    Mv6List = []
    n6List = []
    for i in range(len(MV6[l])):
        x.append(MV6[l][i])
        y.append(n6[l][i])
        Mv6List.append([MV6[l][i]])
        n6List.append([n6[l][i]])
        c.append(tempList1[l])
    conc1.append(np.append(Mv6List, n6List, axis=1))
linesN6 = LineCollection([conc1[l] for l in telmod1],linewidth=1,cmap=cm.winter)
res6 = ax5.scatter(x,y,s=20,c=c,marker='d',cmap=cm.winter, linewidths=0.5, label = '$\mathrm{n_{SF}=6\ cm^{-3}}$')

epsilon1=[0.1,0.3,0.5,0.7,0.9]
arr = np.array([epsilon1[l] for l in telmod1])
linesN6.set_array(arr)
ax5.add_collection(linesN6)

conc2 = []
x = []
y = []
c = []
for k in telmod2:
    Mv50List = []
    n50List = []
    for m in range(len(MV50[k])):
        Mv50List.append([MV50[k][m]])
        x.append(MV50[k][m])
        n50List.append([n50[k][m]])
        y.append(n50[k][m])
        c.append(tempList2[k])
    conc2.append(np.append(Mv50List, n50List, axis=1))
linesN50 = LineCollection([conc2[k] for k in telmod2],linewidth=1,cmap=cm.autumn)
res50 = ax5.scatter(x,y,s=20,c=c,marker='^',cmap=cm.autumn, linewidths=0.5, label='$\mathrm{n_{SF}=50\ cm^{-3}}$')
epsilon2=[0.3,0.5,0.7,0.9]
arr = np.array([epsilon2[k] for k in telmod2])
linesN50.set_array(arr)
ax5.add_collection(linesN50)

## Mb_Z
# Observational data from Grebel
ax4.plot(dsph.getData('M_V'), dsph.getData('[Fe/H]'), '.', c='y', markersize=6, label='$\mathrm{Observations}$', zorder=0)
ax4.plot(de.getData('M_V'), de.getData('[Fe/H]'), '.', c='y', markersize=6, zorder=0)
ax4.plot(dirr.getData('M_V'), dirr.getData('[Fe/H]'), '.', c='y', markersize=6, zorder=0)
ax4.plot(dirrdsph.getData('M_V'), dirrdsph.getData('[Fe/H]'), '.', c='y', markersize=6, zorder=0)

def convert(vals):
    global np
    return 4.83 - 2.5*vals

# Observational data from Sharina
ax4.plot(convert(dsphShar.getData('L_V')), dsphShar.getData('[Fe/H]'), '.', c='y', markersize=6, zorder=0)
ax4.plot(convert(dirrShar.getData('L_V')), dirrShar.getData('[Fe/H]'), '.', c='m', markersize=6, zorder=0)
ax4.plot(convert(dirrdsphShar.getData('L_V')), dirrdsphShar.getData('[Fe/H]'), '.', c='m', markersize=6, zorder=0)

# Observational data from Lianou
ax4.plot(dsphLian.getData('M_V'), dsphLian.getData('[Fe/H]'), '.', c='y', markersize=6, zorder=0)
ax4.plot(dirrdsphLian.getData('M_V'), dirrdsphLian.getData('[Fe/H]'), '.', c='m', markersize=6, zorder=0)

# Simulation data
ax4.plot(n1_e1.getData('M_V'), n1_e1.getData('[Fe/H]'),'h--', color='k', markersize=4, label='$\mathrm{n_{SF}=0.1,\ \epsilon_{FB}=0.1}$', zorder=0 )

conc1 = []
x = []
y = []
c = []
for l in telmod1:
    Mv6List = []
    FeH6List = []
    for i in range(len(MV6[l])):
        x.append(MV6[l][i])
        y.append(FeH6[l][i])
        Mv6List.append([MV6[l][i]])
        FeH6List.append([FeH6[l][i]])
        c.append(tempList1[l])
    conc1.append(np.append(Mv6List, FeH6List, axis=1))
linesN6 = LineCollection([conc1[l] for l in telmod1],linewidth=1,cmap=cm.winter)
res6 = ax4.scatter(x,y,s=20,c=c,marker='d',cmap=cm.winter, linewidths=0.5, label = '$\mathrm{n_{SF}=6\ cm^{-3}}$')

epsilon1=[0.1,0.3,0.5,0.7,0.9]
arr = np.array([epsilon1[l] for l in telmod1])
linesN6.set_array(arr)
ax4.add_collection(linesN6)

conc2 = []
x = []
y = []
c = []
for k in telmod2:
    Mv50List = []
    FeH50List = []
    for m in range(len(MV50[k])):
        Mv50List.append([MV50[k][m]])
        x.append(MV50[k][m])
        FeH50List.append([FeH50[k][m]])
        y.append(FeH50[k][m])
        c.append(tempList2[k])
    conc2.append(np.append(Mv50List, FeH50List, axis=1))
linesN50 = LineCollection([conc2[k] for k in telmod2],linewidth=1,cmap=cm.autumn)
res50 = ax4.scatter(x,y,s=20,c=c,marker='^',cmap=cm.autumn, linewidths=0.5, label='$\mathrm{n_{SF}=50\ cm^{-3}}$')
epsilon2=[0.3,0.5,0.7,0.9]
arr = np.array([epsilon2[k] for k in telmod2])
linesN50.set_array(arr)
ax4.add_collection(linesN50)

ax6.set_xlim(-7,-20)
ax3.set_xlim(-7,-20)
ax1.set_ylim(-1.2,1.5)
ax3.set_ylim(0.7,1.4)
ax6.set_ylim(28,10)
#ax2.set_ylim(0.7,2.1)
ax5.set_ylim(-1,4.5)
ax4.set_ylim(-4., 0.)

#for ax in (ax1, ax2):
#    ax.xaxis.label.set_visible(False)
#    for label in ax.get_xticklabels():
#        label.set_visible(False)

for ax in (ax2, ax4):
    ax.get_xticklabels()[0].set_visible(False)

ax1.set_ylabel('$\mathrm{\log_{10}(R_{e})\ [kpc]}$')
#ax2.set_xlabel('$\mathrm{M_{V}}$')
ax3.set_xlabel('$\mathrm{M_{V}}$')
ax1.set_xlabel('$\mathrm{M_{V}}$')
ax6.set_xlabel('$\mathrm{M_{V}}$')
ax5.set_xlabel('$\mathrm{M_{V}}$')
ax4.set_xlabel('$\mathrm{M_{V}}$')
#ax2.set_ylabel('$\mathrm{\log_{10}(\sigma)\ [km/s]}$')
ax3.set_ylabel('$\mathrm{V-I}$')
ax6.set_ylabel('$\mathrm{\mu_{0}}$')
ax5.set_ylabel('$\mathrm{n}$')
ax4.set_ylabel('$\mathrm{[Fe/H]}$')
 
visual.fig.subplots_adjust(left = 0.08, right = 0.90, bottom = 0.04, top = 0.97,hspace = 0.13, wspace=0.2)

name = "bigFigNew.eps"
figurepath = name
visual.finalize(name=figurepath, dpi=200, show=False)
