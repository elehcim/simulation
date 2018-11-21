# Adapted from https://gist.github.com/sonium0/65cde09ec21d573281d1ce2ca33a6262
# Plot logtime in realtime using bokeh and tail -f
# Tested with python 3.5 and bokeh 0.12.4
# OSX/Linux only
# usage:
# 1. run 'bokeh serve'
# 2. run 'python3.5 main.py logfile.csv'

import sys
import datetime
import asyncio
import asyncio.subprocess

from bokeh.models import ColumnDataSource, HoverTool, SaveTool

from bokeh.client import push_session
from bokeh.plotting import figure, curdoc

import re
from parsers.parse_info import pattern

# def push(step, dt):
#     source.stream(dict(x=[step], y=[dt]), 1000)


def push(d):
    source.stream(d, 500)


@asyncio.coroutine
def update():
    # create = asyncio.create_subprocess_exec('tail', '-f', '-n' , '+1', sys.argv[-1],
    create = asyncio.create_subprocess_exec('tail', '-f', sys.argv[-1],
                                            stdout=asyncio.subprocess.PIPE)
    proc = yield from create
    while True:
        # Read one line of output
        data = yield from proc.stdout.readline()
        line = data.decode('ascii').rstrip()

        # print(line)
        if line:
            d = re.match(pattern, line).groupdict()
            # line = data.decode('ascii').rstrip()
            # line = line.split(', ')
            # line format:
            # label,unix-timestamp,fair_price,spread
            # push(d['step'], d['dt'])
            for k, v in d.items():
                d[k] = [float(v)]
            # print(d)
            push(d)


hover = HoverTool(tooltips=[
    ("Time", "@t"),
    ("dt", "@dt"),
    ("Redshift", "@z")
    ])

log_plot = figure(plot_width=800,
                  plot_height=400,
                  y_axis_type='log',
                  tools="xpan,xwheel_zoom,xbox_zoom,reset", #[hover, SaveTool()],
                  title="Simulation timestep")

log_plot.x_range.follow = "end"

source = ColumnDataSource(data=dict(step=[], dt=[], t=[], z=[]))
# p = figure()
# p.xaxis.formatter=DatetimeTickFormatter()
l = log_plot.line(x='t', y='dt', source=source)

# open a session to keep our local document in sync with server
session = push_session(curdoc(), session_id='main')

session.show(log_plot) # open the document in a browser

loop = asyncio.get_event_loop()
loop.run_until_complete(update())