import os
import glob
import pandas as pd
import numpy as np
import bokeh
from astropy.table import Table
from bokeh.layouts import layout, row, column
from bokeh.plotting import figure, show, save, output_notebook, output_file
from bokeh.models import CustomJS, Slider, ColumnDataSource, Range1d, Plot, Span, Legend, LegendItem, ColorBar, LinearColorMapper
from bokeh.models.glyphs import ImageURL
from bokeh.models.widgets import Panel, Tabs
from bokeh.palettes import Spectral6, Viridis256
from bokeh.transform import linear_cmap

from bokeh.resources import INLINE
# output_notebook(resources=INLINE)


def get_maps_img(img_dir):
    """Produce a dictionary with maps as value and index as keys.
    To be then added to the dataframe with pd.Series(dict).
    a typical map name is: `maps_snapshot_????_*_????.png`
                                                 ^   ^
                                                -8  -4
    """
    print("Found {} maps".format(len(maps)))
    maps = sorted(glob.glob(os.path.join(img_dir, 'maps_img', 'maps_*.png')))
    d = dict()
    for m in maps:
        idx = int(m[-8:-4])-1
        d[idx] = m
    return d



print(bokeh.__version__)

SSAM_DIR = '/home/michele/sim/analysis/ssam/derotated/m69p2_no_vel'

TBL = '/home/michele/sim/analysis/ng_ana/data/mb.69002_p200_a800_r600.fits'

# output_file(os.path.join(SSAM_DIR, 'ssam.html'), mode='inline', title="Specific Stellar Angular Momentum - " + SSAM_DIR)
output_file('ssam.html', mode='inline', title="Specific Stellar Angular Momentum - " + SSAM_DIR)

WINDOW=20
TOOLS='pan,wheel_zoom,reset,hover'

# get_data
ssam_tbl = Table.read(os.path.join(SSAM_DIR, 'maps_data_v_w10_r200_n30.fits'))
ssam_tbl.remove_column('lambda_prof')
ssam_df = ssam_tbl.to_pandas()

df = Table.read(TBL).to_pandas()

# Add maps:
d = get_maps_img(SSAM_DIR)
df['maps'] = pd.Series(d)
# df['maps'] = df['maps'].fillna('')
# TODO if SFR is NaN fix it
# df.sfr = np.ones(len(df))

# Merge data
for col in ('time', 'r_eff_kpc', 'r_eff_kpc3d', 'lambda_r', 'ellip', 'theta', 'n', 'Lx', 'Ly', 'Lz', 'mag_v'):  # TODO compare these two 'r_eff_kpc', 'r_eff_kpc3d'
    df[col] = ssam_df[col]


# Do averages
# df['lambda_r_mean'] = df.lambda_r.rolling(window=WINDOW).agg(np.mean)
# df['ellip_mean'] = df.ellip.rolling(window=WINDOW).agg(np.mean)
for col in ('lambda_r', 'ellip', 'n'):
    df[col+'_mean'] = df[col].rolling(window=WINDOW, win_type='gaussian', center=True).mean(std=5)


print(df.head())

# Visualization part:
source = ColumnDataSource(df)

p = figure(x_range=(0,1), y_range=(0,1), tools='', plot_height=430, plot_width=1150, sizing_mode='stretch_both')
im = p.image_url(url=[source.data['maps'][0]], x=-0.11, y=1, w=None, h=None, anchor="top_left")

p.toolbar.logo = None
p.toolbar_location = None
p.xaxis.visible = None
p.yaxis.visible = None
p.xgrid.grid_line_color = None
p.ygrid.grid_line_color = None
p.outline_line_alpha = 0

COMMON_FIG_ARGS = dict(plot_width=450, plot_height=350, tools=TOOLS)


##### MAPS

s1 = figure(**COMMON_FIG_ARGS)
# s1.line(source.data['time'], source.data['lambda_r_mean'], color="navy")
r1 = s1.multi_line(xs=[source.data['time']]*2, ys=[source.data['lambda_r_mean'], source.data['lambda_r']],
                  color=['navy', 'green'], line_alpha=[1, 0.4])
legend1 = Legend(items=[
    LegendItem(label='lambda_mean ' + str(WINDOW), renderers=[r1], index=0),
    LegendItem(label="lambda_r", renderers=[r1], index=1),
])
legend1.location = "top_left"
legend1.border_line_alpha=0.2
s1.add_layout(legend1)

s1.y_range.start = 0
circle = s1.circle(x=source.data['time'][WINDOW], y=source.data['lambda_r_mean'][WINDOW], color='red')
s1.title.text = "Specific Stellar Angular Momentum"
s1.xaxis.axis_label = 'time'
s1.yaxis.axis_label = 'lambda_mean ' + str(WINDOW)


#### REFF

s2 = figure(**COMMON_FIG_ARGS, x_range=s1.x_range)
r2 = s2.multi_line(xs=[source.data['time']]*2, ys=[source.data['r_eff_kpc'], source.data['r_eff_kpc3d']],
                  color=['navy', 'green'])
# s2.line(x=source.data['time'], y=source.data['r_eff_kpc3d'], color='green')
legend2 = Legend(items=[
    LegendItem(label="r_eff", renderers=[r2], index=0),
    LegendItem(label="3D r_eff", renderers=[r2], index=1),
])
legend2.location = "top_left"
s2.add_layout(legend2)

s2.y_range.start = 0
s2.title.text = "Effective radius"
s2.xaxis.axis_label = 'time'
s2.yaxis.axis_label = 'r_eff'

vline = Span(dimension='height', line_color='red')
s2.add_layout(vline)


#### ORBIT

s3 = figure(**COMMON_FIG_ARGS, match_aspect=True)

pos = s3.line(source.data['x'], source.data['y'])
cross = s3.cross(x=source.data['x'][0], y=source.data['y'][0], color='red', size=12)
s3.cross(x=0, y=0, color='green', size=5)
s3.title.text = "Orbit"
s3.xaxis.axis_label = 'x [kpc]'
s3.yaxis.axis_label = 'y [kpc]'


#### MAG

s4 = figure(**COMMON_FIG_ARGS, x_range=s1.x_range)

pos = s4.line(source.data['time'], source.data['mag_v'])
circle_mag = s4.circle(x=source.data['time'][0], y=source.data['mag_v'][0], color='red')
s4.y_range.flipped = True
s4.title.text = "Magnitude"
s4.xaxis.axis_label = 'time'
s4.yaxis.axis_label = 'Magnitude (V band)'


#### ELLIP

# tab1 = Panel(child=p1, title="xy")
# tab2 = Panel(child=p2, title="size")

# tabs = Tabs(tabs=[ tab1, tab2 ])

s5 = figure(**COMMON_FIG_ARGS, x_range=s1.x_range)

l5a = s5.line(source.data['time'], source.data['ellip_mean'], color='navy')
l5b = s5.line(source.data['time'], source.data['ellip'], line_alpha=0.4, color='green')
circle_ell_t = s5.circle(x=source.data['time'][WINDOW], y=source.data['ellip_mean'][WINDOW], color='red')
s5.y_range.start = 0
s5.title.text = "Ellipticity"
s5.xaxis.axis_label = 'time'
s5.yaxis.axis_label = 'ellip_mean ' + str(WINDOW)
legend5 = Legend(items=[
    LegendItem(label='ellip_mean ' + str(WINDOW), renderers=[l5a], index=0),
    LegendItem(label="ellip", renderers=[l5b], index=1),
])
legend5.location = "top_left"
legend5.border_line_alpha=0.2
s5.add_layout(legend5)


#### ELLIP-LAMBDA
s6 = figure(**COMMON_FIG_ARGS)

y = df.time.values
# mapper = LinearColorMapper(palette=Viridis256, low=min(y), high=max(y))
mapper = linear_cmap(field_name='time', palette=Viridis256, low=min(y), high=max(y))

# s6.circle('ellip_mean', 'lambda_r_mean', fill_color={'field': 'time', 'transform': mapper}, source=source)
el_c = s6.circle('ellip_mean', 'lambda_r_mean', line_color=mapper, fill_color=mapper, fill_alpha=1, size=8, source=source)

circle_ell_l = s6.circle(x=source.data['ellip_mean'][WINDOW], y=source.data['lambda_r_mean'][WINDOW],
                         size=el_c.glyph.size+2, fill_color=None, line_color='red')
color_bar = ColorBar(color_mapper=mapper['transform'], width=12, location=(0, 0))
s6.add_layout(color_bar, 'right')
s6.xaxis.axis_label = 'ellip_mean ' + str(WINDOW)
s6.yaxis.axis_label = 'lambda_mean ' + str(WINDOW)


#### SFR

s7 = figure(**COMMON_FIG_ARGS, x_range=s1.x_range)

s7.line(source.data['time'], source.data['sfr'])
circle_sfh = s7.circle(x=source.data['time'][0], y=source.data['sfr'][0], color='red')
s7.title.text = "SFH"
s7.xaxis.axis_label = 'time'
s7.yaxis.axis_label = 'SFR [Msol/yr]'


##### SERSIC

s8 = figure(**COMMON_FIG_ARGS, x_range=s1.x_range)

l8a = s8.line(source.data['time'], source.data['n_mean'], color='navy')
l8b = s8.line(source.data['time'], source.data['n'], line_alpha=0.4, color='green')
circle_n = s8.circle(x=source.data['time'][WINDOW], y=source.data['n_mean'][WINDOW], color='red')
s8.y_range.start = 0
s8.title.text = "Sersic Index"
s8.xaxis.axis_label = 'time'
s8.yaxis.axis_label = 'n_mean ' + str(WINDOW)
legend8 = Legend(items=[
    LegendItem(label='n_mean ' + str(WINDOW), renderers=[l8a], index=0),
    LegendItem(label="n", renderers=[l8b], index=1),
])
legend8.location = "top_left"
legend8.border_line_alpha=0.2
s8.add_layout(legend8)



#### GENERAL
for _plot in [s1, s2, s3, s4, s5, s6, s7, s8]:
    _plot.toolbar.autohide = True

#### CALLBACK

cb = CustomJS(args=dict(span=vline, im=im, source=source, circle=circle.glyph,
                        cross=cross.glyph, circle_mag=circle_mag.glyph,
                        cet=circle_ell_t.glyph, cel=circle_ell_l.glyph,
                        csfr=circle_sfh.glyph, cn=circle_n.glyph), code="""
    circle.x=source.data['time'][cb_obj.value];
    circle.y=source.data['lambda_r_mean'][cb_obj.value];
    circle_mag.x=source.data['time'][cb_obj.value];
    circle_mag.y=source.data['mag_v'][cb_obj.value];
    cet.x=source.data['time'][cb_obj.value];
    cet.y=source.data['ellip_mean'][cb_obj.value];
    cel.x=source.data['ellip_mean'][cb_obj.value];
    cel.y=source.data['lambda_r_mean'][cb_obj.value];

    csfr.x=source.data['time'][cb_obj.value];
    csfr.y=source.data['sfr'][cb_obj.value];

    cn.x=source.data['time'][cb_obj.value];
    cn.y=source.data['n_mean'][cb_obj.value];

    cross.x=source.data['x'][cb_obj.value];
    cross.y=source.data['y'][cb_obj.value];
    span.location = source.data['time'][cb_obj.value];
    im.data_source.data['url'] = [source.data['maps'][cb_obj.value]];
    im.data_source.change.emit();
""")

# slider.js_on_change('value', cb)
slider = Slider(start=0, end=len(df)-1, value=0, step=1, title="image number", callback=cb)

l = layout([p, slider, row(column(s1, s2), column(s3, s4), column(s5, s6), column(s7, s8))])

save(l)