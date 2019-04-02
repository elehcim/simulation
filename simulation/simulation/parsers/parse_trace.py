import sys
import pandas as pd
import numpy as np

def parse_dens_trace(fname='dens_temp_trace.txt'):
    df = pd.read_csv(fname, delim_whitespace=True, header=1, names=['step', 't', 'vel', 'x','y','z','rho','temp'])
    df['r'] = np.sqrt(df.x**2 + df.y**2 + df.z**2)
    return df

def parse_trace(fname='trace.txt', trace_version=1):
    # All.TracePos[0],All.TracePos[1],All.TracePos[2],
          # All.TraceVel[0],All.TraceVel[1],All.TraceVel[2],
          # All.AdhocAcc[0],All.AdhocAcc[1],All.AdhocAcc[2]);
    print('Parsing a version {} trace file'.format(trace_version))
    if trace_version==1:
        df = pt_v1(fname)
    elif trace_version==2:
        df = pt_v2(fname)

    df['r'] = np.sqrt(df.x**2 + df.y**2 + df.z**2)
    df['v'] = np.sqrt(df.vx**2 + df.vy**2 + df.vz**2)
    df['a'] = np.sqrt(df.ax**2 + df.ay**2 + df.az**2)
    return df

def pt_v1(fname):
    return pd.read_csv(fname, delim_whitespace=True, header=4,
        names=['t', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'ax', 'ay', 'az'])


def pt_v2(fname):
    return pd.read_csv(fname, delim_whitespace=True, header=4,
        names=['t', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'ax_adhoc', 'ay_adhoc', 'az_adhoc',
               'ax', 'ay', 'az', 'omega_x', 'omega_y', 'omega_z', 'phi', 'step', 'dt', 'redshift'])


if __name__ == '__main__':
    parse_trace(sys.argv[1])