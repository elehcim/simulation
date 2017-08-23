import matplotlib.pyplot as plt
import os
import numpy as np
import scipy as sp
import math
import matplotlib
from hyplot.plot import PFigure
from hyplot.visual.PRunData import PRunData, PFileData
import chyplot
import enums

matplotlib.rc('font', size=13)

#Tables for the number of RGB stars per 1000 Msol for a particle with a certain age, both for normal stars and Pop 3 stars in billion years

agesPop2 = [0.00316227766017, 0.00344454795152, 0.00375201416997, 0.00408692534689, 0.00445173126604, 0.00484910038303, 0.00528193934437, 0.00575341424882, 0.00626697380646, 0.00682637456514, 0.00743570838856, 0.00809943238715, 0.0088224015206, 0.0096099041106, 0.0104677005234, 0.0114020653054, 0.0124198330797, 0.0135284485394, 0.0147360209037, 0.0160513832346, 0.017484157048, 0.0190448226929, 0.0207447960122, 0.022596511846, 0.0246135149899, 0.0268105592707, 0.0292037154669, 0.0318104888623, 0.034649947292, 0.0377428606186, 0.041111852658, 0.0447815666663, 0.048778845599, 0.0531329284591, 0.0578756641732, 0.063041744557, 0.0686689580771, 0.0747984662628, 0.0814751047916, 0.0887477114502, 0.0966694833691, 0.105298366145, 0.114697477699, 0.124935569962, 0.136087531783, 0.148234936716, 0.161466639707, 0.175879427049, 0.191578724344, 0.208679367662, 0.227306443536, 0.247596203936, 0.269697062916, 0.29377068222, 0.319993153797, 0.348556287861, 0.379669015932, 0.413558919117, 0.450473892796, 0.490683959915, 0.53448324613, 0.582192131251, 0.634159592737, 0.690765758369, 0.752424686782, 0.819587396191, 0.892745163463, 0.972433117677, 1.05923415444, 1.15378319963, 1.25677185367, 1.36895345043, 1.49114856763, 1.62425102916, 1.76923444318, 1.92715932373, 2.09918084817, 2.28655730694, 2.49065930764, 2.71297980065, 2.95514499962, 3.21892627681, 3.50625312019, 3.81922724712, 4.16013797782, 4.53147898113, 4.93596651503, 5.37655929531, 5.85648013777, 6.37923953224, 6.94866132087, 7.56891066844, 8.24452452946, 8.98044483473, 9.78205464019, 10.6552175025, 11.6063203696, 12.6423202988, 13.7707953468, 15.0]

nRGBsPop2 = [2.0265928267e-96, 0.0415643666089, 0.0476819494392, 0.0499839649844, 0.0539812249985, 0.0600816254442, 0.0659567826436, 3.14933998568e-14, 0.014594040174, 0.00595944943027, 0.00704017423502, 0.00629418492754, 0.00735265442973, 0.00702737749785, 0.0090025984562, 0.0118910651631, 0.0167093352359, 0.0336303108195, 0.0534799383644, 0.0752545629205, 0.107134717069, 0.124643994902, 0.162230153562, 0.17345023982, 0.19359274554, 0.202002623681, 0.22033528447, 0.217574329707, 0.245735812097, 0.233116562941, 0.262594004275, 0.297232378027, 0.288382749239, 0.328576298637, 0.346523572021, 0.36566608496, 0.387330768998, 0.424971770584, 0.461284830491, 0.497645245606, 0.533800762743, 0.564940649908, 0.598808894368, 0.632065716926, 0.681957201853, 0.739497866213, 0.770581110501, 0.796653196945, 0.838603135092, 0.882588116554, 0.889228403125, 0.935067622311, 0.921045263135, 0.916466050286, 0.916463099707, 0.891446262125, 0.898514722113, 0.884282317901, 0.866645610609, 0.892079716842, 0.858258522864, 0.847267693996, 0.819831044166, 0.845585767005, 0.858740522302, 0.871830389294, 0.912688150659, 0.961635391434, 1.03019356391, 1.0428008681, 0.981094308499, 0.834821473013, 1.37369607651, 1.6357086316, 1.78098553986, 1.91772620422, 2.08093748326, 2.22186895342, 2.66614049712, 3.52610484618, 4.45050197047, 4.04539757268, 3.79146879058, 3.6662755411, 3.24615691361, 3.15402818544, 3.03825961548, 3.0684759681, 3.07179983459, 3.09647085574, 3.05811295047, 3.10676407968, 3.0188818674, 3.08166227424, 2.95363368129, 2.99524219601, 2.89788240339, 2.98829597498, 3.08434890072, 3.03037862315]

agesPop3 = [7.0, 7.16326530612, 7.32653061224, 7.48979591837, 7.65306122449, 7.81632653061, 7.97959183673, 8.14285714286, 8.30612244898, 8.4693877551, 8.63265306122, 8.79591836735, 8.95918367347, 9.12244897959, 9.28571428571, 9.44897959184, 9.61224489796, 9.77551020408, 9.9387755102, 10.1020408163, 10.2653061224, 10.4285714286, 10.5918367347, 10.7551020408, 10.9183673469, 11.0816326531, 11.2448979592, 11.4081632653, 11.5714285714, 11.7346938776, 11.8979591837, 12.0612244898, 12.2244897959, 12.387755102, 12.5510204082, 12.7142857143, 12.8775510204, 13.0408163265, 13.2040816327, 13.3673469388, 13.5306122449, 13.693877551, 13.8571428571, 14.0204081633, 14.1836734694, 14.3469387755, 14.5102040816, 14.6734693878, 14.8367346939, 15.0]

nRGBsPop3 = [0.000259933756365, 0.000254655670083, 0.000249035236802, 0.00024309098949, 0.000236842041345, 0.000230308085006, 0.000223509389661, 0.000216466796649, 0.000209201712622, 0.0002017361004, 0.000194092467134, 0.00018629384947, 0.00014430945021, 0.000102298548442, 0.000100025807825, 9.76974321841e-05, 9.53174209724e-05, 9.28898198446e-05, 9.04187188451e-05, 8.79082504211e-05, 8.53625872492e-05, 8.27859398616e-05, 8.01825540601e-05, 7.75567081032e-05, 7.49127096535e-05, 7.2254892467e-05, 6.95876128096e-05, 6.69152455824e-05, 6.42421801358e-05, 6.15728157528e-05, 5.89115567791e-05, 5.62628073757e-05, 5.36309658701e-05, 5.10204186769e-05, 4.84355337594e-05, 4.58806536006e-05, 4.33600876482e-05, 4.08781041971e-05, 3.84389216689e-05, 3.60466992448e-05, 3.3705526803e-05, 2.78693920928e-05, 2.01438930385e-05, 1.98636420234e-05, 2.05650190295e-05, 2.11151421373e-05, 2.15220047706e-05, 2.17935714688e-05, 2.19377719886e-05, 2.19624951148e-05]

#File to where the calculated properties are written to

RGBfile=open("results/RGBsRe.dat", "w")

RGBfile.writelines('#Sim\tnRGBpop2\tnRGBpop3\t[<Fe/H>]\t<[Fe/H]>')

# Read in a list of simulations and some of their previously calculated properties. 
# If you do not have these yet, you can use some dummy values just to test the script
# props = PFileData("results/generalproperties.dat")
# simulations = props.getData('Sim')
# R30s = props.getData('R_30')
# MVs = props.getData('M_V')

simulations = [60003]
R30s = [1]
MVs = [-15]

#Function to plot the observational luminosity-metallicity relation
def plot_LV_FeH(ax, xMin=3, xMax=9, yMax=-0.5, yMin=-2.5):
	global PFileData, np

	obsdir = "/home/michele/sim/analysis/observationalData/"

	KI13_MWdSph = PFileData(os.path.join(obsdir, "13Kirby_MWdSph.dat"))
	KI13_M31dSph = PFileData(os.path.join(obsdir, "13Kirby_M31dSph.dat"))
	KI13_LGdIrr = PFileData(os.path.join(obsdir, "13Kirby_LGdIrr.dat"))

	# Observational data from Kirby 2013 (spectroscopic metallicities)
	ax.plot((KI13_LGdIrr.getData('log(L_V/L_sun)')), KI13_LGdIrr.getData('[Fe/H]'), ".",  markersize=6, zorder=0, color='0.5')
	ax.plot((KI13_M31dSph.getData('log(L_V/L_sun)')), KI13_M31dSph.getData('[Fe/H]'), ".",  markersize=6, zorder=0, color='0.5')
	ax.plot((KI13_MWdSph.getData('log(L_V/L_sun)')), KI13_MWdSph.getData('[Fe/H]'), ".",  markersize=6, zorder=0, color='0.5')

	ax.set_xlim(xMin, xMax)
	ax.set_ylim(yMin, yMax)

	#Fit to the total data (Kirby et al. 2013)
	slope = 0.21
	intercept = -2.84
	fitrange = np.arange(xMin,xMax+0.1,0.1)

	fit = np.array([intercept+f*slope for f in fitrange])
	ax.plot(fitrange, fit, '--', color='0.5', zorder=0)

	ax.set_xlabel('$\log_{10}(L_V/\mathrm{L_\odot})$')
	ax.set_ylabel('$\mathrm{\langle[Fe/H]\\rangle}_\mathrm{RGB}$')


def w_average(prop, weighing_factor, log=False):
	#Calculate the weighted average of a certain property.
	#The first argument is a list of the properties for each stellar particle (e.g. the [Fe/H]s)
	#2nd argument is the property to which it is weighted (e.g. number of RGB stars in a stellar particle)
	#3rd optional argument determines the difference between the log of the average or the average of the log (<[Fe/H]> vs [<Fe/H>]). 
	#Kirby 2014 uses the average of the log, so standard is False

	fehSol = -2.756
	global math
	tot = sum(weighing_factor)
	if not log:
		weighed_prop = [p*w for p, w in zip(prop, weighing_factor)]
		return sum(weighed_prop)/tot
	else:
		weighed_prop = [10**(p+fehSol)*w for p, w in zip(prop, weighing_factor)]		
		return math.log10(sum(weighed_prop)/tot)-fehSol

fig = plt.figure(FigureClass = PFigure.PFigure, figsize=(6, 3))
ax1 = fig.add_my_subplot(111)

plot_LV_FeH(ax1) #Plot the observations

for sim, R30, MV in zip(simulations, R30s, MVs):

	#Read the data
	sim = int(sim)
	dr = chyplot.CDataGadget(sim)
	fdir = "/home/michele/sim/sim{}".format(sim) #Adjust to your own path with simulations

	dr.setPrefix( fdir )
	dr.checkFilesPresent() # set the first and last dump
	dr.set_file( 100) #) dr.lastDump())
	print dr.lastDump()
	data = dr.readFile()

	x,y,z = data.rcom(True, enums.T_star, 0, 0, 0, True)
	data.vcom(True, enums.T_star)
	data.convertUnits()

	#Only look at stars within a certain area (R30 is the radius where the surface density drops below 30 mag/"^2
	visitorR = chyplot.cglobals.plmap.getSecond("radius")
	data.applyLimits(visitorR, 0, R30, enums.T_all) 

	#Divide the data into Pop3 and Pop2 stellar particles
	visitor = chyplot.cglobals.plmap.getSecond("[Fe/H]")
	dataPop3 = data.limitsCopy(visitor, -99, -5, enums.T_star)
	dataPop2 = data.limitsCopy(visitor, -5, 100, enums.T_star)

	#Get lists of the ages and masses of the stars
	agesPop2Data = dataPop2.getDataArray(enums.T_star, chyplot.cglobals.plmap.getSecond('birthtime'), False)
	agesPop2Data = [data.time() - a for a in agesPop2Data]
	massesPop2Data = dataPop2.getDataArray(enums.T_star, chyplot.cglobals.plmap.getSecond('initialMass'), False)

	agesPop3Data = dataPop3.getDataArray(enums.T_star, chyplot.cglobals.plmap.getSecond('birthtime'), False)
	agesPop3Data = [data.time() - a for a in agesPop3Data]
	massesPop3Data = dataPop3.getDataArray(enums.T_star, chyplot.cglobals.plmap.getSecond('initialMass'), False)

	#Get the number of RGB stars per stellar particle, using the lists at the beginning of the script
	RGBsPop2 = np.interp(agesPop2Data, agesPop2, nRGBsPop2)
	RGBsPop3 = np.interp(agesPop3Data, agesPop3, nRGBsPop3)

	RGBsPop2 = [r*m*1e3 for r, m in zip(RGBsPop2, massesPop2Data)]
	RGBsPop3 = [r*m*1e3 for r, m in zip(RGBsPop3, massesPop3Data)]

	#Since the number of RGB stars in a Pop3 stellar particle is neglible, we only use the metallicity of the Pop2s.
	metalPop2 = dataPop2.getDataArray(enums.T_star, chyplot.cglobals.plmap.getSecond('[Fe/H]'), False)
	
	weightedFeH = w_average(metalPop2, RGBsPop2)

	galaxy = chyplot.CGalaxy(dataPop2)
	Z_lum, FeH_lum, MgFe_lum = galaxy.getTotalStellarMetallicityLum()

	RGBfile.writelines('\n{}\t{}\t{}\t{}\t{}'.format(sim, sum(RGBsPop2), sum(RGBsPop3), np.log(weightedFeH), weightedFeH))

	LV = -0.4*(MV-4.8)

	print LV, weightedFeH

	ax1.plot(LV, weightedFeH, 'b', marker='o', markersize=6, markeredgecolor='w', markeredgewidth=0.2)



directory = '/home/michele/sim/results/metallicity/' #Adjust
name = 'metallicityRGB.pdf'

	 
fig.subplots_adjust(left = 0.1, right = 0.975, bottom = 0.15, top = 0.95,wspace = 0.2)

# plt.show()
# fig.finalize(name=directory + name, dpi=300, show=False)













