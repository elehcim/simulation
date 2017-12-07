#!/bin/env python3
import re
import pandas as pd
from itertools import zip_longest


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

patterns = [r'Step= (?P<step>{})  t= (?P<t>{})  dt= (?P<dt>{})'.format(i, g, g),
            r'Nf= (?P<Nf>{})  total-Nf= (?P<tot_Nf>{})  ex-frac= (?P<ex>{})  iter= (?P<iter>{})'.format(i, i, g, i),
            r'work-load balance: (?P<wl>{})  max=(?P<max>{}) avg=(?P<avg>{}) PE0=(?P<PE0>{})'.format(g, g, g, g),
            r'particle-load balance: (?P<plb>{})'.format(g),
            r'max. nodes: (?P<max_nodes>{}), filled: (?P<filled>{})'.format(i, g),
            r'part/sec=(?P<part_sec>{}) \| (?P<boh>{})  ia/part=(?P<ia_part>{})'.format(g, g, g)
]

def parse_timings(fname="timings.txt"):
    with open("timings.txt") as f:
        content = f.read().splitlines()
    frames = list()
    counter = 0
    for lines in grouper(content, 7):
        counter += 1
        line_dict = {}
        # print(lines)
        for i, line in enumerate(lines[:-1]):
            # print(line)
            d = re.match(patterns[i], line).groupdict()
            # print(d)
            line_dict.update(d)
        if counter % 1000 == 0:
            print(line_dict)
        frames.append(pd.DataFrame([line_dict]))

    return pd.concat(frames, ignore_index=True)

# Line0
# m = re.match(r'Step= (\d+)  t= (\d+(\.\d*)?|\.\d+)', "Step= 0  t= 0.9")
# m = re.match(r'Step= (?P<step>\d+)  t= (?P<t>\d+(\.\d*)?|\.\d+)', line[0])