#=============================================================================
# A python script to extract and save the properties ()  of star particles
# from a simulation's snapshot.
# run as : '>>> python starprops.py <sim_no> <snapshot_no>'
#=============================================================================

import numpy as np
import enums
import os
import chyplot
import sys

# sys.argv[0] is ths name of the script itself
# simulation = sys.argv[1]
# snapshot = int( sys.argv[2] )

simulation = sys.argv[1]
snapshot = int( sys.argv[2] )

onlyRe = True
Re = 1

## fdir ="/run/user/1000/gvfs/smb-share:server=files.ugent.be,share=srathi/shares/vakgroep/astro_dwarfgalaxies/"

# fdir = "/run/user/1000/gvfs/smb-share\:server\=files.ugent.be\,share\=srathi/shares/vakgroep/astro_dwarfgalaxies/"
fdir = "/home/michele/sim/"
dr = chyplot.CDataGadget()
# dr.setPrefix( fdir + "rpverbek/sim{}/sim{}".format( simulation[:-1]+"0", simulation ) ) #for hidden simulations (ending in 0)

## dr.setPrefix( fdir + "rpverbek/sim{}".format( simulation ) )
dr.setPrefix( fdir + "sim{}".format( simulation ) )
dr.checkFilesPresent()
dr.set_file( snapshot )

data = dr.readFile()
data.rcom( True, enums.T_star, 0, 0, 0, True )
data.rotate( enums.T_gas, 2, True )

# Choosing particles within one half light radius
if onlyRe:
	visitor = chyplot.cglobals.plmap.getSecond("dist_to_z")
	data.applyLimits(visitor, 0, Re, enums.T_star)

# Throwing out Pop III star particles
visitor2 = chyplot.cglobals.plmap.getSecond("[Fe/H]")
data.applyLimits(visitor2, -5, 100, enums.T_star)

data.convertUnits()
massesStars = data.getDataArray(enums.T_star, chyplot.cglobals.plmap.getSecond('initialMass'), True)
birthtimes = data.getDataArray(enums.T_star, chyplot.cglobals.plmap.getSecond('birthtime'), True)
ages = [data.time() - b for b in birthtimes]
metals = data.getDataArray(enums.T_star, chyplot.cglobals.plmap.getSecond('metallicity'), True)
alphas = data.getDataArray(enums.T_star, chyplot.cglobals.plmap.getSecond('[Mg/Fe]'), True)
iron = data.getDataArray(enums.T_star, chyplot.cglobals.plmap.getSecond('[Fe/H]'), True)

massesStars = np.array( massesStars )*1e6
print massesStars.sum(), 'Msol'
ages = np.array( ages )
metals = np.array( metals )
alphas = np.array( alphas )
iron = np.array( iron )

exit(0)

starprops = np.vstack((massesStars, ages, metals, alphas, iron)).T

header = '{}	{}	{}	{}	{} \n{}	{}	{}	{}	{}'.format('massPart', 'Age', 'metal', '[Mg/Fe]', '[Fe/H]', 'Msol', 'Gyrs', 'no_unit', 'dex', 'dex')

# fdir = fdir + 'srathi/sim{}/'.format( simulation )
# if not os.path.exists( fdir ):
# 	os.system('mkdir -p {}'.format( fdir ) )
# f = fdir + 'starprops_{}.dat'.format( simulation )
# np.savetxt(f, starprops, delimiter = '	', header = header )
