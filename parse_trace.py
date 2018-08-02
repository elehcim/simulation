import sys
import pandas as pd
import numpy as np

def parse_dens_trace(fname='dens_temp_trace.txt'):
    df = pd.read_csv(fname, delim_whitespace=True, header=1, names=['step', 't', 'vel', 'x','y','z','rho','temp'])
    df['r'] = np.sqrt(df.x**2 + df.y**2 + df.z**2)
    return df

def parse_trace(fname='trace.txt'):
    # All.TracePos[0],All.TracePos[1],All.TracePos[2],
          # All.TraceVel[0],All.TraceVel[1],All.TraceVel[2],
          # All.AdhocAcc[0],All.AdhocAcc[1],All.AdhocAcc[2]);
    df = pd.read_csv(fname, delim_whitespace=True, header=4, names=['t', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'ax', 'ay', 'az'])
    df['r'] = np.sqrt(df.x**2 + df.y**2 + df.z**2)
    return df


if __name__ == '__main__':
    parse_trace(sys.argv[1])