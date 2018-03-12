import chyplot
import enums

simulation = 10000 #simulation number
snapshot = 66 #snapshot number (e.g. 1)


#read in the snapshot
dr = chyplot.CDataGadget(simulation)
fdir = "/home/michele/sim/sim{}".format(simulation)

dr.setPrefix( fdir )
dr.checkFilesPresent() # set the first and last dump
dr.set_file( snapshot )
data = dr.readFile()

data.rcom(True, enums.T_star, 0, 0, 50, True)
data.vcom(True, enums.T_star)
data.convertUnits()

#get the total mass of all star particles
particleProperty = chyplot.cglobals.plmap.getSecond('mass')
massesStars = data.getDataArray(enums.T_star, particleProperty, True)
totalStarMass = sum(massesStars)*1e6

print "Total stellar mass is {:.0f} Msol".format(totalStarMass)
