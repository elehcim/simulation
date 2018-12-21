#!/usr/bin/env python3
import re
import pandas as pd
from itertools import zip_longest

# From Gadget2 Documentation
# 7.2 Performance statistics of the gravitational tree
#
# In the file specified by TimingsFile, you get some statistics about the performance of the
# gravitational force computation by the tree algorithm, which is usually (but not always) the
# main sink of computational time in a simulation. A typical output may look like this:
#
#   Step= 482 t= 1.29199 dt= 0.00292969
#   Nf= 10992 total-Nf= 0016981370
#   work-load balance: 1.10997 max=0.618782 avg=0.557475 PE0=0.618782
#   particle-load balance: 1.02133
#   max. nodes: 15202, filled: 0.422278
#   part/sec=9858.74 | 8881.96 ia/part=1004.71 (0)
#   Step= 483 t= 1.29492 dt= 0.00292969
#   Nf= 60000 total-Nf= 0017041370
#   work-load balance: 1.01677 max=2.51368 avg=2.47221 PE0=2.43074
#   particle-load balance: 1.02133
#   max. nodes: 15202, filled: 0.422278
#   part/sec=12134.9 | 11934.7 ia/part=805.744 (0)
#
#
# The first line of the block generated for each step informs you about the number of the current
# timestep, the current simulation time, and the system timestep itself (i.e. the time difference to
# the last force computation).
# `Nf` gives the number of particles that receive a force computation in the current timestep,
# while `total-Nf` gives the total number of force computations since the simulation was started.
# `work-load balance` gives the work-load balance in the actual tree walk, i.e. the largest
# required time of all the processors divided by the average time.
# The number given by `max` is the maximum time (in seconds) among the processors spent for the tree walk,
# while `avg` gives the average, and `PE0` the time for processor 0.
# `particle-load balance` gives the maximum number of particles per processor divided by the
# average number, i.e. it measures the memory imbalance. An upper bound for this number is PartAllocFactor.
# However, for technical reasons, the domain decomposition always leaves room for a few par-
# ticles, hence this upper bound is never exactly reached. `max.nodes` informs about the maxi-
# mum number of internal tree nodes that are used among the group of processors, and the number
# behind `filled` gives this quantity normalised to the total number of allocated tree nodes, it
# hence gives the fraction to which the tree storage is filled. Finally, `part/sec` measures the
# raw force speed in terms of tree-force computations per processor per second. The first num-
# ber gives the speed that would be achieved for perfect work-load balance. Due to work-load
# imbalance, the actually achieved average effective force speed will always be lower in practice,
# and this number is given after the vertical line.
# `ia/part` gives the average number of particle node interactions required to compute the force
# for each of the active particles.
# The number in parenthesis is only non-zero if the Ewald correction is used to obtain periodic boundaries.
# It then gives the average number of Ewald correction-forces evaluated per particle.



# from here: https://stackoverflow.com/a/434328/1611927
def chunker(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))


# from here https://stackoverflow.com/a/434411/1611927
def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


regex_float = r'(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?'
regex_int = r'\d+'
g = regex_float
i = regex_int
ordered_columns = ["step", "t", "dt",
                   "Nf", "tot_Nf", "ex", "iter",
                   "wl", "max", "avg", "PE0",
                   "plb",
                   "max_nodes", "filled",
                   "ps", "ps_real", "ia_part"]

patterns = [r'Step= (?P<step>{})  t= (?P<t>{})  dt= (?P<dt>{})'.format(i, g, g),
            r'Nf= (?P<Nf>{})  total-Nf= (?P<tot_Nf>{})  ex-frac= (?P<ex>{})  iter= (?P<iter>{})'.format(i, i, g, i),
            r'work-load balance: (?P<wl>{})  max=(?P<max>{}) avg=(?P<avg>{}) PE0=(?P<PE0>{})'.format(g, g, g, g),
            r'particle-load balance: (?P<plb>{})'.format(g),
            r'max. nodes: (?P<max_nodes>{}), filled: (?P<filled>{})'.format(i, g),
            r'part/sec=(?P<ps>{}) \| (?P<ps_real>{})  ia/part=(?P<ia_part>{})'.format(g, g, g)
]


def parse_timings(fname="timings.txt"):
    print("Reading...")
    with open(fname) as f:
        content = f.read().splitlines()
    print("Parsing...")
    frames = list()
    for lines in grouper(content, 7):
        line_dict = {}
        for i, line in enumerate(lines[:-1]):
            d = re.match(patterns[i], line).groupdict()
            line_dict.update(d)
        frames.append(line_dict)
    print('Converting to DataFrame...')
    df = pd.DataFrame(frames)
    # reorder columns
    df = df[ordered_columns]
    print('Casting types...')
    # Convert to numerics
    for col in df.columns:
        df[col] = pd.to_numeric(df[col])
    return df

if __name__ == '__main__':
    df = parse_timings()
    print(df.head())