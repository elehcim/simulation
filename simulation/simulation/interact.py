import numpy as np
import matplotlib.pylab as plt
import pynbody
from ipywidgets import (interactive, HBox, VBox, Checkbox, Layout,
                        IntSlider, FloatRangeSlider, ToggleButtons, SelectMultiple)

_RHO_MAX = 2e-1
_RHO_MIN = 5e-4


def _available_keys():  # FIXME do it in an automamtic way. But seems that self.profiles is not useful and better not use it.
    # keys = set(self[0].g.loadable_keys()).union(set(sim.profiles[0]['g'].derivable_keys()))
    # return keys
    keys = ['E_circ', 'acce_norm', 'vel_norm', 'Q', 'X', 'beta', 'density', 'density_enc', 'dyntime', 'fesp', 'fourier',
    'g_spherical', 'j_circ', 'j_phi', 'j_theta', 'jtot', 'kappa', 'magnitudes', 'mass', 'mass_enc',
    'mgsp', 'omega', 'pattern_frequency', 'pot', 'p', 'rho', 'rotation_curve_spherical', 'sb',
    'smooth', 'temp', 'u', 'v_circ', 'vel', 'zsph', 'vr', 'vr_disp']

    return keys

# FIXME I need to find a unified way to interact...


def interact_gas(sim, rho_min=_RHO_MIN, rho_max=_RHO_MAX, step=1e-5):

    def k(i, velocity_proj, sfh, cog, vrange, width, resolution):
        sim.plot_gas(i, velocity_proj=velocity_proj, sfh=sfh, cog=cog,
                                vmin=vrange[0], vmax=vrange[1],
                                width=width, resolution=resolution)

    _vminmax = FloatRangeSlider(
        value=[rho_min, rho_max],
        min=rho_min,
        max=rho_max,
        step=step,
        description='Rho:',
        continuous_update=False,
        orientation='horizontal',
        readout=True,
        readout_format='1.0e',
    )
    _snap_slider = IntSlider(min=0,max=len(sim)-1,step=1,value=0, continuous_update=False, description='Snap:')
    _width_slider = IntSlider(min=5,max=1000,step=10,value=20, continuous_update=False, description='Width (kpc):')
    _res_slider = IntSlider(min=100,max=1000,step=100,value=200, continuous_update=False, description='Resol. (pix):')
    _proj = Checkbox(value=False, description='Velocity projection')
    _sfh = Checkbox(value=True, description='SFH')
    _traj = Checkbox(value=True, description='COG traj.')

    w = interactive(k,
                    i=_snap_slider,
                    velocity_proj=_proj,
                    sfh=_sfh,
                    cog=_traj,
                    vrange=_vminmax,
                    width=_width_slider,
                    resolution=_res_slider)
    control_column_size = 4
    b = VBox([HBox([VBox(w.children[0:control_column_size]), VBox(w.children[control_column_size:-1])],
             layout=Layout(display='flex', width='150%')), w.children[-1]])
    return b

def interact(sim, rho_min=_RHO_MIN, rho_max=_RHO_MAX, step=1e-5):
    def k(i, velocity_proj, sfh, cog, vrange, width, starsize, resolution):
        sim.plot_gas_and_stars(i, velocity_proj=velocity_proj, sfh=sfh, cog=cog,
                                vmin=vrange[0], vmax=vrange[1], starsize=starsize,
                                width=width, resolution=resolution)

    from ipywidgets import interactive, HBox, VBox, Layout, FloatRangeSlider, IntSlider, Checkbox, FloatSlider

    _vminmax = FloatRangeSlider(
        value=[rho_min, rho_max],
        min=rho_min,
        max=rho_max,
        step=step,
        description='Rho:',
        continuous_update=False,
        orientation='horizontal',
        readout=True,
        readout_format='1.0e',
    )
    _snap_slider = IntSlider(min=0,max=len(sim)-1,step=1,value=0, continuous_update=False, description='Snap:')
    _width_slider = IntSlider(min=5,max=1000,step=10,value=20, continuous_update=False, description='Width (kpc):')
    _res_slider = IntSlider(min=100,max=1000,step=100,value=200, continuous_update=False, description='Resol. (pix):')
    _proj = Checkbox(value=False, description='Velocity projection')
    _sfh = Checkbox(value=True, description='SFH')
    _traj = Checkbox(value=True, description='COG traj.')
    _starsize = FloatSlider(min=0.1, max=1000, value=1, continuous_update=False, description='Starsize. (kpc):')

    w = interactive(k,
                    i=_snap_slider,
                    velocity_proj=_proj,
                    sfh=_sfh,
                    cog=_traj,
                    vrange=_vminmax,
                    width=_width_slider,
                    starsize=_starsize,
                    resolution=_res_slider)
    b = VBox([HBox([VBox(w.children[0:4]), VBox(w.children[4:8])],
             layout=Layout(display='flex', width='150%')), w.children[-1]])
    return b


def interact_profiles(sim, default='u', eps=0.03, keys=_available_keys(), add_keys=None, selection=True,
                          _snap_slider=None, offset=0, log=False, **kwargs):
    # if keys is None:
    #     keys = _available_keys()
    if add_keys is not None:
        keys += add_keys

    if _snap_slider is None:
        _snap_slider = IntSlider(min=0,max=len(sim)-1,step=1,value=0, continuous_update=False, description='Snap:')
    _family_choice = ToggleButtons(options=['g','s'], value='g')

    def create_vminmax(kwargs):
        min = kwargs['min'] if 'min' in kwargs else 0.0
        max = kwargs['max'] if 'max' in kwargs else 1000.0

        _vminmax = FloatRangeSlider(
            value=[min, max],
            min=min,
            max=max,
            step=10,
            description='radius:',
            continuous_update=False,
            orientation='horizontal',
            readout=True,
            readout_format='.0f',
        )
        return _vminmax

    _vminmax = create_vminmax(kwargs)

    def create_varoptions(default=default, selection_type=SelectMultiple):
        return selection_type(options=keys, value=default)

    def k(i, family, y, vrange=None):
        if vrange is not None:
            kwargs.update({'min': vrange[0], 'max':vrange[1]})

        i = offset + i
        from pynbody import units
        snap = sim[i]
        snap['eps'] = pynbody.array.SimArray(eps*np.ones_like(sim[i]['mass']), units.kpc)
        pynbody.analysis.halo.center(snap)

        f = getattr(snap, family)
        # Add norm of vector quantities
        if 'acce' in f.loadable_keys():
            f['acce_norm'] = np.linalg.norm(f['acce'], axis=1)
        if 'vel' in f.loadable_keys():
            f['vel_norm'] = np.linalg.norm(f['vel'], axis=1)

        p = pynbody.analysis.profile.Profile(f, **kwargs)

        # if ax is None:
        fig, ax = plt.subplots(1, figsize=(6,4))
    #     print(sim.profiles[i])
    #     print(p)
    #     if not isinstance(y, (list, tuple)):
    #         y = tuple(y)
    #     for _y in y:

        ax.plot(p['rbins'], p[y])
        if log:
            ax.semilogy()
        snap_time_gyr = snap.properties['time'].in_units("Gyr")
        ax.set_xlabel("r ({})".format(p['rbins'].units))
        try:
            ax.set_ylabel("{} ({})".format(y, r"${}$".format(getattr(p[y],'units','').latex())))
        except AttributeError:
            ax.set_ylabel("{} ({})".format(y, getattr(p[y],'units','')))
        title = '{}   ($t={:5.2f}$ Gyr, snap={})'.format(y, snap_time_gyr, i)
        ax.set_title(title)

    w = interactive(k,
                    i=_snap_slider, family=_family_choice, vrange=_vminmax,
                    y=default if not selection else create_varoptions(default, selection_type=ToggleButtons))
    return w