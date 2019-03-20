import sys
import pandas as pd

def parse_sfc(fname='sfc.dat'):
    df = pd.read_csv(fname, delim_whitespace=True, header=0, names=['t', 'ndens', 'ndiv','nentr','nphyscond','nrandom','nbelowJeans','nfarbelowJeans','TotN_gas'])
    return df

# int ndens = 0,      // How many particles (per node per timestep) satisfy the density threshold criterion
#     ndiv = 0,       // How many particles satisfy the negative divergence criterion
#     nentr = 0,      // How many particles satisfy the temperature threshold criterion
#     nphyscond = 0,  // How many particles satisfy all the condition above
#     nrandom = 0,    // How many particles satisfy the stochastic star formation
#     nbelowJeans = 0,     // How many particles are below the Jeans temperature
#     nfarbelowJeans = 0;  // How many particles are below 0.9 * Jeans temperature
# 8.386      31  130166   86688      31       0       0       0 130166




if __name__ == '__main__':
    parse_sfc(sys.argv[1])