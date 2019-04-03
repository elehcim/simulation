import sys
import pandas as pd
import numpy as np

COLUMNS_FORMAT = {
1 : ['t', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'ax_adhoc', 'ay_adhoc', 'az_adhoc'],
2 : ['t', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'ax_adhoc', 'ay_adhoc', 'az_adhoc',
     'ax', 'ay', 'az', 'omega_x', 'omega_y', 'omega_z', 'phi', 'step', 'dt', 'redshift']
}

HEADER_LINES = {1: 3, 2: 4}

def parse_dens_trace(fname='dens_temp_trace.txt'):
    df = pd.read_csv(fname, delim_whitespace=True, header=1, names=['step', 't', 'vel', 'x','y','z','rho','temp'])
    df['r'] = np.sqrt(df.x**2 + df.y**2 + df.z**2)
    return df


def get_trace_version(fname):
    with open(fname) as fp:
        for i, line in enumerate(fp):
            if i == 4:
                break
    if 'TraceAcc' in line:
        return 2
    else:
        return 1


def parse_trace(fname='trace.txt', force_trace_version=None):
    if force_trace_version is not None:
        trace_version = force_trace_version
    else:
        trace_version = get_trace_version(fname)
    print('Parsing a version {} trace file'.format(trace_version))
    df = pd.read_csv(fname, delim_whitespace=True,
                            header=HEADER_LINES[trace_version],
                            names=COLUMNS_FORMAT[trace_version])

    df['r'] = np.sqrt(df.x**2 + df.y**2 + df.z**2)
    df['v'] = np.sqrt(df.vx**2 + df.vy**2 + df.vz**2)
    df['a_adhoc'] = np.sqrt(df.ax_adhoc**2 + df.ay_adhoc**2 + df.az_adhoc**2)
    return df


if __name__ == '__main__':
    parse_trace(sys.argv[1])