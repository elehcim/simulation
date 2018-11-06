import os
import glob
import pandas as pd
import numpy as np
import bokeh
from bokeh.layouts import layout, row, column
from bokeh.plotting import figure, show, save, output_notebook, output_file
from bokeh.models import CustomJS, Slider, ColumnDataSource, Range1d, Plot, Span, Legend, LegendItem
from bokeh.models.glyphs import ImageURL
from bokeh.models.widgets import Panel, Tabs

from bokeh.resources import INLINE
# output_notebook(resources=INLINE)
output_file('ssam.html', mode='inline')

print(bokeh.__version__)

DATA_DIR = 'ssam.69p2.xy.sph.new'
# DATA_DIR_S = 'ssam.69p2.side.sph.new'
SIM_PATH = '/home/michele/sim/MySimulations/hi_osc/mb.69002_p200_a800_r600/out'

def get_data(filename, generate=False, save_data=True):
    if not generate:
        return pd.read_pickle(filename)
    maps = glob.glob(os.path.join(DATA_DIR, 'maps_*.png'))
    maps.sort()

    df = pd.read_csv(os.path.join(DATA_DIR, "tl.dat"), sep=' ', skiprows=1,
        names=['time', 'lambda_r', 'ellip', 'theta', 'r_eff_kpc', 'r_eff_kpc3d', 'Lx', 'Ly', 'Lz', 'mag_v'])
    df['maps'] = maps
    print(df.head())


    from parse_trace import parse_trace
    sim_path = '/home/michele/sim/MySimulations/hi_osc/mb.69002_p200_a800_r600/out'
    # sim = simulation.Simulation(sim_path)
    trace = parse_trace(os.path.join(sim_path, 'trace.txt'))

    # Sync trace time and snapshot time
    # Right is true to keep the last index otherwise it will result in NaN
    locations = np.digitize(df.time, trace.t, right=True)
    print(locations.shape)

    for f in ['r', 'v', 'a', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'ax', 'ay', 'az']:
        df[f] = trace[f].loc[locations].values

    print(df.vx.head())

    ## Do rolling averages
    df['lambda_r_mean'] = df.lambda_r.rolling(window=20).agg(np.mean)
    df['ellip_mean'] = df.ellip.rolling(window=20).agg(np.mean)

    if save_data:
        df.to_pickle(filename)
    return df

df = get_data('d1.pkl', generate=False, save_data=True)

# Visualization part:
source = ColumnDataSource(df)


s1 = figure(plot_width=350, plot_height=250)
s1.line(source.data['time'], source.data['lambda_r_mean'], color="navy")

vline = Span(dimension='height', line_color='red')

p = figure(x_range=(0,1), y_range=(0,1), tools='', plot_height=430, plot_width=1150, sizing_mode='stretch_both')
im = p.image_url(url=[source.data['maps'][0]], x=-0.01, y=1, w=None,h=None, anchor="top_left")

vline = Span(dimension='height', line_color='red')

circle = s1.circle(x=source.data['time'][0], y=source.data['lambda_r'][0])
# im.data_source.data['url'] = [source.data['maps'][5]]


p.toolbar.logo = None
p.toolbar_location = None
p.xaxis.visible = None
p.yaxis.visible = None
p.xgrid.grid_line_color = None
p.ygrid.grid_line_color = None
p.outline_line_alpha = 0

############

s2 = figure(x_range=s1.x_range, plot_width=350, plot_height=250)
# r = s2.multi_line(xs=[source.data['time']]*2, ys=[source.data['r_eff_kpc'], source.data['r_eff_kpc3d']], color=['blue', 'green'])
r = s2.multi_line(xs=[source.data['time']]*2, ys=[source.data['r_eff_kpc'], source.data['r_eff_kpc3d']], color=['blue', 'green'])
# s2.line(x=source.data['time'], y=source.data['r_eff_kpc3d'], color='green')

##########

s3 = figure(plot_width=350, plot_height=250, match_aspect=True)

pos = s3.line(source.data['x'], source.data['y'])
cross = s3.cross(x=source.data['x'][0], y=source.data['y'][0], color='red')
s3.add_layout(cross)
######################

s4 = figure(plot_width=350, plot_height=250)

pos = s4.line(source.data['time'], source.data['mag_v'])
circle_mag = s4.circle(x=source.data['time'][0], y=source.data['mag_v'][0])
s4.add_layout(circle_mag)
######################

# tab1 = Panel(child=p1, title="xy")
# tab2 = Panel(child=p2, title="size")

# tabs = Tabs(tabs=[ tab1, tab2 ])


# p = figure()
# r = p.multi_line([[1,2,3], [1,2,3]], [[1,3,2], [3,4,3]],
#                  color=["orange", "red"], line_width=4)

legend = Legend(items=[
    LegendItem(label="r_eff", renderers=[r], index=0),
    LegendItem(label="3D r_eff", renderers=[r], index=1),
])
legend.location = "top_left"
s2.add_layout(legend)


cb = CustomJS(args=dict(span=vline, im=im, source=source, circle=circle.glyph,
                        cross=cross.glyph, circle_mag=circle_mag.glyph), code="""
    circle.x=source.data['time'][cb_obj.value];
    circle.y=source.data['lambda_r_mean'][cb_obj.value];
    circle_mag.x=source.data['time'][cb_obj.value];
    circle_mag.y=source.data['mag_v'][cb_obj.value];
    cross.x=source.data['x'][cb_obj.value];
    cross.y=source.data['y'][cb_obj.value];
    span.location = source.data['time'][cb_obj.value];
    im.data_source.data['url'] = [source.data['maps'][cb_obj.value]];
    im.data_source.change.emit();
""")

# slider.js_on_change('value', cb)
slider = Slider(start=0, end=len(df.maps)-1, value=0, step=1, title="image number", callback=cb)
s1.add_layout(vline)

s1.add_layout(circle)

l = layout([p, slider, row(column(s1, s2), column(s3, s4))])

save(l)
