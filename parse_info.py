import sys
import re
import pandas as pd
from itertools import zip_longest

# Example of a line:
# Begin Step 35508, Time: 0.202861, Redshift: 9.49992, Systemstep: 1.59878e-06

regex_float = r'(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?'
regex_int = r'\d+'
g = regex_float
i = regex_int
ordered_columns = ["step", "t", "z", "dt"]

pattern = r'Begin Step (?P<step>{}), Time: (?P<t>{}), Redshift: (?P<z>{}), Systemstep: (?P<dt>{})'.format(i, g, g, g)
# pc = re.compile(pattern)

def parse_info(fname="info.txt"):
    with open(fname) as f:
        content = f.read().splitlines()
    frames = list()
    for line in content:
        if line:
            line_dict = {}
            d = re.match(pattern, line).groupdict()
            line_dict.update(d)
            frames.append(line_dict)
    df = pd.DataFrame(frames)
    # reorder columns
    df = df[ordered_columns]
    # Convert to numerics
    for col in df.columns:
        df[col] = pd.to_numeric(df[col])
    return df

if __name__ == '__main__':
    df = parse_info(sys.argv[-1])