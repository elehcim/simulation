import sys
import pandas as pd
from itertools import islice


# header from Gadget2/run.c:343
header = ["tot", "grav", "hydro", "domain",
"potential", "drift", "kick",
"write_snap", "tree_walk", "tree_build",
"comm", "imbalance", "sph",
"comm_sph", "imb_sph", "neigh_adj",
"PM", "peano", "add_phys", "star"]

header_new = ["tot", "grav", "hydro", "domain",
"potential", "drift", "kick",
"write_snap", "tree_walk", "tree_build",
"comm", "imbalance", "sph",
"comm_sph", "imb_sph", "neigh_adj",
"PM", "peano", "add_phys", "star_hydro", "star_fb", "star"]

def parse_cpu(fname="cpu.txt", new_format=False):
    with open(fname, "r") as f:
        nums = list()
        # Skip odd lines, keep only even
        # usage of islice from here: https://stackoverflow.com/a/16022845/1611927
        for line in islice(f, 1, None, 2):
            nums.append([float(n) for n in line.split()])
    h = header_new if new_format else header
    df = pd.DataFrame(nums, columns=h)
    return df

if __name__ == '__main__':
    df = parse_cpu(sys.argv[-1])
    print(df.mean())
    print(df.describe())