# from bokeh.plotting import figure, show, output_notebook
# import bokeh
# print(bokeh.__version__)
# p = figure(x_range=(0,1), y_range=(0,1))

# url = '/home/michele/sim/analysis/ssam.69p2.xy.sph/maps_snapshot_0002_v_w10_r400_a20.png'
# # url = 'http://users.ugent.be/~mmastrop/maps_snapshot_0004_v_w10_r400_a20.png'
# p.image_url(url=[url], x=0, y=0, w=1, h=0.3, anchor="bottom_left")
# # show(p)
# show(p)
import os
import matplotlib.pyplot as plt
# import simulation
import numpy as np
import pandas as pd
from parse_trace import parse_trace


from bokeh.layouts import layout

from bokeh.models import CustomJS, ColumnDataSource, Slider
from bokeh.plotting import figure, output_file, show

output_file('image.html')

source = ColumnDataSource(data=dict(url=['maps_snapshot_0002_v_w10_r400_a20.png']))

p = figure(x_range=(0,1), y_range=(0,1))

callbackimg = CustomJS(args=dict(source=source), code="""
    var data = source.data;
    var f = cb_obj.value;
    var old = data['url'][0];
    var padded = "_".concat(("0000" + f).slice(-4)).concat("_");
    var patt1 = /_\d\d\d\d_/;
    var result = old.match(patt1);
    var out = old.replace(result, padded);
    console.log(out);
    data['url'][0] = out;
    source.change.emit();
""")
    #old.replace("1",f.toString(10))

p.image_url('url',source=source, x=0, y=1,w=1,h=1)
slider = Slider(start=2, end=500, value=2, step=2, title="image number", callback=callbackimg)


def linked_panning():
	source = ColumnDataSource(data=dict(url=['maps_snapshot_0002_v_w10_r400_a20.png']))
	a = np.loadtxt('ssam.69p2.xy.sph/tl.dat')

    N = 100
    x = np.linspace(0, 4 * np.pi, N)
    y1 = np.sin(x)
    y2 = np.cos(x)
    y3 = np.sin(x) + np.cos(x)

    s1 = figure(tools=tools)
    s1.circle(x, y1, color="navy", size=8, alpha=0.5)
    s2 = figure(tools=tools, x_range=s1.x_range, y_range=s1.y_range)
    s2.circle(x, y2, color="firebrick", size=8, alpha=0.5)
    s3 = figure(tools='pan, box_select', x_range=s1.x_range)
    s3.circle(x, y3, color="olive", size=8, alpha=0.5)

    vline = Span(location=0, dimension='height', line_color='red', line_width=3)
    s3.add_layout(vline)
    return [s1, s2, s3]

l = layout([slider, p])
show(l)