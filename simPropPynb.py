################################################################################
# Using pynbody, make a table of simulation properties:
# sim_no | M* | r_200 | r_vir | M_vir
# NOT OPTIMIZED
# Author: Shivangee Rathi, Michele Mastropietro
################################################################################
import numpy as np
import sys, glob, os
import pynbody
import pandas as pd

def get_sims( direc ):
    dir_names = list(glob.glob(direc+'/M??sim?????'))
    dir_names.sort()
    print(dir_names)
    sims = []
    for f in dir_names:
        sims.append( int( f[-5:] ) )
    sims = np.array( sims )
    return dir_names, sims

def getM200( r200, h ): # r200 is in kpc
    Msol = 1.988547e30 # solar mass in kg
    kpc = 3.085677581e19 # parsec in m

    H0 = 100 * h # H0 = 100 h km s-1 Mpc-1
    G = 6.67408*10**(-11) # in m3 kg-1 s-2
    return (100*(H0**2)*(r200**3)/G)*(1000**(-2)*kpc/(Msol*10**(-6))) # in Msol

def getSimProps( snap_path ):
    snap = pynbody.load(snap_path)
    snap.properties['omegaM0'] = 0.28
    pynbody.analysis.halo.center( snap ) # FIXME center on the stars? RE-center

    m_star = np.float( snap.s['massform'].in_units('Msol').sum() ) #total stellar mass
    r_eff = np.float( pynbody.analysis.luminosity.half_light_r(snap, band ='v' ).in_units('kpc') ) # effective radius - 3D,; for 2D set 'cylindrical=True'
    r200 = np.float( pynbody.analysis.halo.virial_radius(snap, overden=200, rho_def='critical').in_units('kpc') ) # virial radius -- using overdensity of 200
    M200 = getM200( r200, snap.properties['h'] ) # mass enclosed within r200 in Msol

    return m_star, r_eff, r200, M200

if __name__ == '__main__':

    direc = '/home/michele/sim/MoRIA/M1-10_Verbeke2017/'
    # direc = sys.argv[1]

    dir_names, sims = get_sims( direc )

    simProps = {}

    for sim_dir, sim_name in zip(dir_names, sims):
        snap_path = list(glob.glob(sim_dir+'/snapshot_????'))
        snap_path.sort()
        simProps[sim_name] = getSimProps( snap_path[-1] )

# simProps = np.zeros([len(files), 4])
# for i, sim_num in enumerate(sims):
#     # print(sims[i])
#     # if not os.path.exists( direc + "rpverbek/sim{}".format( sims[i] ) ):
#     #     sim_path = direc + "rpverbek/sim{}0/sim{}/".format( files[i][-5:-1], sims[i] )
#     # else:
#     #     sim_path = direc + "rpverbek/sim{}/".format( sims[i] )
#     sim_path = os.path.join(direc, files[i])
#     snaps = os.listdir( sim_path )
#     snaps.sort()
#     snap_path = sim_path + snaps[-1]

#     # print sims[i], snap
#     simProps[i, :] = getSimProps( snap_path )

    df = pd.DataFrame.from_dict(data=simProps, orient='index')
    columns = ['Mstar_tot', 'r_eff', 'r200', 'M200']
    df.columns = columns

# simPropArr = np.array( simProps )

# sim_names = [np.int(sim) for sim in simPropArr[:, 0]]
# simPropArr[:, 0] = sim_names

# header = '{}	{}	{}	{}	{} \n{}	{}	{}	{}	{}'.format('sim_no', 'Mstar_tot', 'r_eff', 'r200', 'M200', 'no_unit', 'Msol', 'kpc', 'kpc', 'Msol')

    out_file = 'simProp_IAU.dat'
# if os.path.exists( out_file ):
# 	os.system('rm -p {}'.format( out_file ) )
    header = ['Mstar_tot (Msol)', 'r_eff (kpc)', 'r200 (kpc)', 'M200 (Msol)']
# df.to_csv(out_file, header=header, sep='\t')
    df.to_latex(out_file, header=header)
# np.savetxt(out_file, simPropArr, delimiter = '	', header = header )
