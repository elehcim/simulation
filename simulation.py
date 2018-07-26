import os
import sys
import numpy as np
import matplotlib.pylab as plt
import csv
import logging
import pynbody
from snap_io import load_moria_sim_and_kicked, load_moria, load_kicked, load_sim
from util import np_printoptions
from analyze_sumfiles import get_sumfile
from multiprocessing import Pool, Process, Queue
from functools import lru_cache


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
logger.addHandler(logging.StreamHandler())

def mass_resolution(snap):
    return (snap['mass'].sum()/len(snap)).in_units("Msol")

def velocity_projection(snap):
    v_x, v_y, v_z = snap['vel'].mean(axis=0)
    v_xy = np.linalg.norm([v_x, v_y])
    alpha = np.sign(v_y) * np.arccos(v_x/v_xy) * 180.0/np.pi
    theta = np.arctan(v_z/v_xy) * 180.0/np.pi
    return alpha, theta

def _cog_queue(snap, q):
    mass = snap['mass']
    pos = snap['pos']
    tot_mass = mass.sum()
    q.put(np.sum(mass * pos.transpose(), axis=1) / tot_mass)


def my_cog(snap):
    # snap = self.snap_list[i]
    mass = snap['mass']
    pos = snap['pos']
    tot_mass = mass.sum()
    return np.sum(mass * pos.transpose(), axis=1) / tot_mass

def plot_cog(cog, ax_cog=None, cur_snap=None, **kwargs):
    if ax_cog is None:
        fig, ax_cog = plt.subplots(1, figsize=(8,8))
    ax_cog.set_xlabel("x [kpc]")
    ax_cog.set_ylabel("y [kpc]")
    ax_cog.scatter(*cog[:2], **kwargs)
    # Plot current position and center
    if cur_snap is not None:
        ax_cog.scatter(*cog[:2, cur_snap], color="red")
    ax_cog.scatter(0, 0, marker='+', color="b")
    ax_cog.set_title("COG trajectory")
    ax_cog.axis('equal')
    return ax_cog
# def _cog(sim, i):
#     snap = sim.snap_list[i]
#     mass = snap['mass']
#     pos = snap['pos']
#     tot_mass = mass.sum()
#     return np.sum(mass * pos.transpose(), axis=1) / tot_mass

def get_param_used(path):
    d = {}
    path = os.path.expanduser(path)
    if os.path.isdir(path):
        path = os.path.join(path, 'parameters-usedvalues')
    with open(path) as csvfile:
        spamreader = csv.reader(csvfile, delimiter=' ', skipinitialspace=True)
        for line in spamreader:
            try:
                k, v = line
                d[k] = v
            except ValueError as e:
                print("{}. But continuing".format(e), file=sys.stderr)
    return d

class Simulation(object):
    """docstring for Simulation"""
    _times = None
    cog = None
    _rho_max=2e-1
    _rho_min=5e-4;
    _computed_cog = False

    def __init__(self, sim_dir, sim_id=None):
        # TODO join this __init__ with the Moria or Kicked class
        if sim_id is None:
            sim_id = sim_dir
        self.sim_id = sim_id
        self._sim_dir = sim_dir
        self.snap_list = self._load(sim_dir)
        if len(self.snap_list) == 0:
            raise RuntimeError("No snaphots found in {}".format(sim_dir))

        self._centered = np.zeros(len(self.snap_list), dtype=bool)
        self.params = get_param_used(sim_dir) 

    def _load(self, sim_id):
        logger.info("loading simulation: {}".format(sim_id))
        return load_sim(sim_id)

    def _cog(self, i):
        snap = self.snap_list[i]
        mass = snap['mass']
        pos = snap['pos']
        tot_mass = mass.sum()
        return np.sum(mass * pos.transpose(), axis=1) / tot_mass

    def snap(self, idx):
        return self.snap_list[idx]

    @property
    def mass_resolution_gas(self):
        return mass_resolution(self[0].g)

    @property
    def mass_resolution_dm(self):
        return mass_resolution(self[0].dm)

    @property
    def mass_resolution(self):
        return mass_resolution(self[0])
    
    @property
    def mass_resolution_baryons(self):
        snap = self[0]
        return ((snap.g['mass'].sum() + snap.s['mass'].sum())/len(snap)).in_units("Msol")

    @property
    @lru_cache(1)
    def times(self):
        return np.array([snap.properties['time'].in_units('Gyr') for snap in self.snap_list])

    def get_times(self):
        self._times = np.zeros(len(self))
        for i, snap in enumerate(self.snap_list):
            self._times[i] = snap.properties['time'].in_units('Gyr')
        return self._times

    def mass_resolution(self):
        return mass_resolution(self.snap_list[0])

    def __getitem__(self, idx):
        return self.snap_list[idx]

    def __len__(self):
        return len(self.snap_list)

    def __repr__(self):
        return "{}: ({}) {}".format(self.sim_id, len(self), self.snap(0).__repr__())

    @property
    def t_range(self):
        if self._times is None:
            self.get_times()
        return self._times.min(), self._times.max()

    # def create_profiles(self, **kwargs):
    #     '''Create profiles objects'''
    #     self.profiles = list() #[None] * len(self)
    #     for i, snap in enumerate(self.snap_list):
    #         self.profiles.append({})
    #         for f in snap.families():
    #             self.profiles[i].update({f.aliases[-1]:
    #                 pynbody.analysis.profile.Profile(getattr(snap, f.name), **kwargs)})


    # def plot_property(a, k, prop, unit=None, ax=None):
    #     if ax is None:
    #         fig, ax = plt.subplots()
    #     unit = k[prop].unit if unit is None else unit
    #     time = a['time'].to(u.Gyr)
    #     ax.plot(time, a[prop].to(unit), "r", label="MoRIA")
    # #     ax2 = ax.twinx()
    # #     ax2.s
    #     k2 = np.interp(a['time'], k['time'], k[prop].to(unit), left=np.nan)
    #     ax.plot(time[:-1], k2[:-1], label="Kicked")
    #     ax.set_xlabel('time (Gyr)')
    #     ax.set_ylabel("{0} ({1:latex})".format(prop,unit))
    #     ax.set_title(prop)
    #     plt.legend(loc=0)


    def compute_cog(self, use_process=False, use_multiprocess=False,
                    save_cache=False, cache_file=None, cache_dir=".cog/", force=False, verbose=True, family=None):
        """
        Compute the center of gravity of a simulation

        Args:
            save_cache: Save the resultsin a npz file
            cache_file: Name of the cache npz file
            verbose: Print which file is currently read
            family: The subfamily attribute to take into account, if `None` consider the whole snapshot.

        Returns:
            cog: is a 3 column array with the coordinates of the center of gravity positions for each snapshot
        """
        if self._computed_cog and not force:
            logger.info("Center of gravity already computed")
            return

        if cache_file is None:
            cache_file = os.path.join(cache_dir, 
                os.path.basename(self.sim_id) + (".cog.npz" if family is None else ".{}.cog.npz".format(family)))
            if save_cache:
                os.makedirs(cache_dir, exist_ok=True)

        if not force and os.path.isfile(cache_file):
            logger.info("Loading precomputed center of gravity for all the snapshots ({})".format(cache_file))
            self.cog = np.load(cache_file)['cog']
            self._computed_cog = True
            return

        logger.info("Computing center of gravity for all the snapshots")

        self.cog = np.zeros((3, len(self)), dtype=float)

        # if use_multiprocess:  # Not working for now, it is stuck it seems because of a thread lock in SimSnap. 
        #     import multiprocess
        #     pool = multiprocess.Pool(processes=8)
        #     multiple_results = pool.map(my_cog, self.snap_list)
        #     return multiple_results
        # else:
        for i, snap in enumerate(self.snap_list):
            if family is not None:
                snap = getattr(snap, family)

            if verbose:
                logger.info("{:03d} Analysing {} (time {:.4f} Gyr)".format(i, snap.filename, snap.properties['time'].in_units('Gyr')))


            if use_process:
                q = Queue()
                p = Process(target=_cog_queue, args=(snap, q))
                p.start()
                this_cog = q.get()
                p.join()
            else:
                mass = snap['mass']
                pos = snap['pos']
                tot_mass = mass.sum()
                this_cog = np.sum(mass * pos.transpose(), axis=1) / tot_mass

            self.cog[:,i] = this_cog

        if save_cache:
            logger.info("Saving cog file: {}".format(cache_file))
            np.savez(cache_file, cog=self.cog)

        self._computed_cog = True
        return self.cog

    def center(self, i):
        if not self._centered[i]:
            pynbody.analysis.halo.center(self[i])
            self._centered[i] = True

    def _center_all(self):
        for snap in self.snap_list:
            pynbody.analysis.halo.center(snap)

    def plot_cog(self, ax_cog=None, cur_snap=None):
        return plot_cog(self.cog, ax_cog, cur_snap)

    def plot_sfh(self, ax_sfh, snap_time_gyr=None, last_snap=-1, **kwargs):
        if ax_sfh is None:
            fig, ax_sfh = plt.subplots(1, figsize=(8,6))
        # ignore AccuracyWarning that is issued when an integral is zero
        import warnings
        from scipy.integrate.quadrature import AccuracyWarning
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=AccuracyWarning)
            pynbody.plot.stars.sfh(self.snap_list[last_snap], subplot=ax_sfh, **kwargs)
        if snap_time_gyr is not None:
            ax_sfh.axvline(x=snap_time_gyr, linestyle="--")
        # ax_sfh.set_title("SFH")
        ax_sfh.set_xlabel("Time [Gyr]")
        ax_sfh.set_ylabel("SFR [M$_\odot$ yr$^{-1}$]")
        return ax_sfh

    def plot_gas_and_stars(self, i, velocity_proj=False, sfh=False, cog=False, starsize=None, **kwargs):
        """Create figure with gas and star rendering from pynbody"""
        snap = self.snap_list[i]
        snap.g['smooth'] /= 2
        snap.s['smooth'] /= 2
        snap_time_gyr = snap.properties['time'].in_units("Gyr")
        pynbody.analysis.halo.center(snap)

        if velocity_proj:
            """x-axis is aligned with the overall mean velocity of
            the snaphot and the vertical axis is the z axis rotated by the elevation angle
            of the velocity"""
            backup = snap['pos'].copy()
            backup_v = snap['vel'].copy()
            alpha, theta = velocity_projection(snap)
            r1 = snap.rotate_z(alpha)
            r2 = snap.rotate_y(theta)

        fig, (ax_g, ax_s) = plt.subplots(nrows=1, ncols=2, figsize=(16,8))

        # If not provided use a default value for width
        width = kwargs.get("width", 20)
        kwargs.pop("width", None)

        try:
            snap_num = int(snap.filename[-4:])
    #         with np_printoptions(precision=2):
    #             title = '$t={:5.2f}$ Gyr, snap={}\nv = {}'.format(snap.properties['time'].in_units("Gyr"), snap_num, velocity)
            im = pynbody.plot.sph.image(snap.g, qty="rho", units="g cm^-2", subplot=ax_g, #title=title,
                           ret_im=True, show_cbar=False, width=width, **kwargs)
            rgbim = pynbody.plot.stars.render(snap, starsize=starsize, axes=ax_s, width=width, clear=False, plot=False, ret_im=True)
            ax_s.imshow(rgbim[::-1, :], extent=(-width / 2, width / 2, -width / 2, width / 2))
            ax_s.set_xlabel('x [' + str(snap.s['x'].units) + ']')
            ax_s.set_ylabel('y [' + str(snap.s['y'].units) + ']')
            ax_g.set_xlabel('x [' + str(snap.g['x'].units) + ']')
            ax_g.set_ylabel('y [' + str(snap.g['y'].units) + ']')

            fig.tight_layout()  # only plots above are affected
            fig.subplots_adjust(top=0.92, bottom=0.15)
            cbar_ax = fig.add_axes([0.12, 0.07, 0.3, 0.02])
            fig.colorbar(im, cax=cbar_ax, orientation='horizontal').set_label("rho [g cm^-2]")
            if sfh:
                # TODO fix negative position of axes
                #  [left, bottom, width, height]
                ax_sfh = fig.add_axes([0.1,  -0.3, 0.35, 0.26])
                self.plot_sfh(ax_sfh, snap_time_gyr)

            if cog:
                ax_cog = fig.add_axes([0.63,  -0.3, 0.26, 0.26])
                self.plot_cog(ax_cog, i)

            title = '$t={:5.2f}$ Gyr, snap={}'.format(snap_time_gyr, snap_num)
            if velocity_proj:
                with np_printoptions(precision=2):
                    mean_velocity = snap['vel'].mean(axis=0)
                    title+="\nv = {} km/s".format(mean_velocity)
            fig.suptitle(title)
        except Exception as e:
            raise(e)
        finally:
            snap.g['smooth'] *= 2
            snap.s['smooth'] *= 2
            if velocity_proj:
                # revert is costly (7% for each transformation w.r.t. the sph.image function)
                # and does not work when the transformation has been applied on a particle family
                # r2.revert()
                # r1.revert()
                snap['pos'] = backup
                snap['vel'] = backup_v
        return fig

    def plot_gas(self, i, velocity_proj=False, sfh=False, cog=False, **kwargs):
        """Create figure with gas and star rendering from pynbody"""
        snap = self.snap_list[i]
        snap.g['smooth'] /= 2
        snap_time_gyr = snap.properties['time'].in_units("Gyr")
        pynbody.analysis.halo.center(snap)

        if velocity_proj:
            """x-axis is aligned with the overall mean velocity of
            the snaphot and the vertical axis is the z axis rotated by the elevation angle
            of the velocity"""
            backup = snap['pos'].copy()
            backup_v = snap['vel'].copy()
            alpha, theta = velocity_projection(snap)
            r1 = snap.rotate_z(alpha)
            r2 = snap.rotate_y(theta)

        fig, ax_g = plt.subplots(nrows=1, ncols=1, figsize=(16,8))

        # If not provided use a default value for width
        width = kwargs.get("width", 20)
        kwargs.pop("width", None)

        try:
            snap_num = int(snap.filename[-4:])
    #         with np_printoptions(precision=2):
    #             title = '$t={:5.2f}$ Gyr, snap={}\nv = {}'.format(snap.properties['time'].in_units("Gyr"), snap_num, velocity)
            im = pynbody.plot.sph.image(snap.g, qty="rho", units="g cm^-2", subplot=ax_g, #title=title,
                           ret_im=True, show_cbar=False, width=width, **kwargs)
            ax_g.set_xlabel('x [' + str(snap.g['x'].units) + ']')
            ax_g.set_ylabel('y [' + str(snap.g['y'].units) + ']')

            fig.tight_layout() # only plots above are affected
            fig.subplots_adjust(top=0.92, bottom=0.15)
            cbar_ax = fig.add_axes([0.35, 0.07, 0.3, 0.02])
            fig.colorbar(im, cax=cbar_ax, orientation='horizontal').set_label("rho [g cm^-2]")
            if sfh:
                # TODO fix negative position of axes
                #  [left, bottom, width, height]
                ax_sfh = fig.add_axes([0.15,  -0.3, 0.35, 0.26])
                self.plot_sfh(ax_sfh, snap_time_gyr)
            if cog:
                ax_cog = fig.add_axes([0.63,  -0.3, 0.26, 0.26])
                self.plot_cog(ax_cog, i)

            title = '$t={:5.2f}$ Gyr, snap={}'.format(snap_time_gyr, snap_num)
            if velocity_proj:
                with np_printoptions(precision=2):
                    mean_velocity = snap['vel'].mean(axis=0)
                    title+="\nv = {} km/s".format(mean_velocity)
            fig.suptitle(title)
        except Exception as e:
            raise(e)
        finally:
            snap.g['smooth'] *= 2
            if velocity_proj:
                # revert is costly (7% for each transformation w.r.t. the sph.image function)
                # and does not work when the transformation has been applied on a particle family
                # r2.revert()
                # r1.revert()
                snap['pos'] = backup
                snap['vel'] = backup_v
        return fig

    # FIXME I need to find a unified way to interact...
    def interact_gas(self, rho_min=None, rho_max=None, step=1e-5):
        if rho_min is None:
            rho_min = self._rho_min
        if rho_max is None:
            rho_max = self._rho_max

        def k(i, velocity_proj, sfh, cog, vrange, width, resolution):
            self.plot_gas(i, velocity_proj=velocity_proj, sfh=sfh, cog=cog,
                                    vmin=vrange[0], vmax=vrange[1],
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
        _snap_slider = IntSlider(min=0,max=len(self)-1,step=1,value=0, continuous_update=False, description='Snap:')
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

    def interact(self, rho_min=None, rho_max=None, step=1e-5):
        if rho_min is None:
            rho_min = self._rho_min
        if rho_max is None:
            rho_max = self._rho_max

        def k(i, velocity_proj, sfh, cog, vrange, width, starsize, resolution):
            self.plot_gas_and_stars(i, velocity_proj=velocity_proj, sfh=sfh, cog=cog,
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
        _snap_slider = IntSlider(min=0,max=len(self)-1,step=1,value=0, continuous_update=False, description='Snap:')
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

    def _available_keys():  # FIXME do it in an automamtic way. But seems that self.profiles is not useful and better not use it.
        # keys = set(self[0].g.loadable_keys()).union(set(sim.profiles[0]['g'].derivable_keys()))
        # return keys
        keys = ['E_circ', 'acce_norm', 'vel_norm', 'Q', 'X', 'beta', 'density', 'density_enc', 'dyntime', 'fesp', 'fourier', 
        'g_spherical', 'j_circ', 'j_phi', 'j_theta', 'jtot', 'kappa', 'magnitudes', 'mass', 'mass_enc',
        'mgsp', 'omega', 'pattern_frequency', 'pot', 'p', 'rho', 'rotation_curve_spherical', 'sb', 
        'smooth', 'temp', 'u', 'v_circ', 'vel', 'zsph', 'vr', 'vr_disp']

        return keys

    def interact_profiles(self, default='u', eps=0.03, keys=_available_keys(), add_keys=None, selection=True,
                          _snap_slider=None, offset=0, log=False, **kwargs):
        # if keys is None:
        #     keys = _available_keys()
        if add_keys is not None:
            keys += add_keys
        from ipywidgets import interactive, IntSlider, ToggleButtons, SelectMultiple, FloatRangeSlider, Select

        if _snap_slider is None:
            _snap_slider = IntSlider(min=0,max=len(self)-1,step=1,value=0, continuous_update=False, description='Snap:')
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
            snap = self[i]
            snap['eps'] = pynbody.array.SimArray(eps*np.ones_like(self[i]['mass']), units.kpc)
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
        #     print(self.profiles[i])
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

    def plot_stars(self, snap):
        """Wrapper around pynbody.plot.stars.render using correct smoothing length"""
        pass

def plot_gas(sim, i, faceon=False, **kwargs):
    snap = sim.snap_list[i]
    snap.g['smooth'] /= 2
    pynbody.analysis.halo.center(snap)
    if faceon:
        pynbody.analysis.angmom.faceon(snap)
    try:
        img = pynbody.plot.sph.image(snap.g, qty="rho", units="g cm^-2", **kwargs)
    except Exception as e:
        raise(e)
    finally:
        snap.g['smooth'] *= 2
    return img


class MoriaSim(Simulation):
    """docstring for MoriaSim"""
    _sf_tidal_folder = "/home/michele/sim/MySimulations/Moria8Gyr_tidal/results/sumfiles/"
    _sf_moria = "/home/michele/sim/MoRIA/results/sumfiles/"

    def __init__(self, sim_id, kicked=False):
        self.sim_id = sim_id
        self.kicked = kicked
        # self.snap_list = load_kicked(sim_id) if kicked else load_moria(sim_id)
        # super(MoriaSim, self).__init__(sim_id)
        self._load(sim_id, kicked)
        sumfile_path = os.path.join(self._sf_moria, sim_id + ".dat")
        if os.path.isfile(sumfile_path):
            logger.info("Getting sumfile: {}".format(sumfile_path))
            self.sumfile = get_sumfile(os.path.join(self._sf_moria, sim_id + ".dat"))
        else:
            logger.info("No sumfile found")
        self._centered = np.zeros(len(self.snap_list), dtype=bool)

    def _load(self, sim_id, kicked=False):
        logger.info("loading simulation: {}".format(sim_id))
        self.snap_list = load_kicked(sim_id) if kicked else load_moria(sim_id)
        # Remove boxsize which complicates the plotting
        for i, snap in enumerate(self.snap_list):
            if i==0:
                self.boxsize = snap.properties.pop('boxsize', None).copy()
            snap.properties.pop('boxsize', None)

    def clear_snap_list(self):
        import gc
        # import resource
        # mem_1 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        # print(mem_1)
        del self.snap_list
        nobs = gc.collect()
        # mem_2 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss  
        # print(mem_2)
        # logger.info("freed {} objects - {:.2f} MB".format(nobs, round((mem_1-mem_2)/1024.0,1)))
        logger.info("freed {} objects".format(nobs))
        self._load(self.sim_id, self.kicked)

    def plot_orbit(self):
        fig, ax = plt.subplots()
        ax.plot(self.sumfile['rcom_x'], self.sumfile['rcom_y'], label="{}".format(self.sim_id))
        ax.set_aspect('equal')
        plt.legend(loc=0)


class BhSim(Simulation):
    """docstring for MoriaSim"""
    _sim_dir = "/home/michele/sim/MySimulations/bh"
    _sf_bh_folder = os.path.join(_sim_dir, "results/sumfiles/")

    def __init__(self, sim_id):
        self.sim_id = sim_id
        # self.snap_list = load_kicked(sim_id) if kicked else load_moria(sim_id)
        # super(MoriaSim, self).__init__(sim_id)
        self._load(sim_id)

    def _load(self, sim_id):
        logger.info("loading simulation: {}".format(sim_id))
        self.snap_list = load_sim(os.path.join(self._sim_dir, sim_id, "out"))
        # Remove boxsize which complicates the plotting
        for i, snap in enumerate(self.snap_list):
            if i==0:
                self.boxsize = snap.properties.pop('boxsize', None).copy()
            snap.properties.pop('boxsize', None)

def time_range_kicked_moria():
    trange = (min(np.min(times_moria), np.min(times_kicked)), max(np.max(times_moria), np.max(times_kicked)))
    trange = (0, min(np.max(times_moria), np.max(times_kicked)))
    bins = np.linspace(*trange, bins_sfr)
    trange, times_moria[-1], times_kicked[-1]

if __name__ == '__main__':
    SIMNUMBER = "69002_p200.0_a600.0_r600.0_c8.15"
    kicked=True
    sim = MoriaSim(SIMNUMBER, kicked)
    s = sim.snap_list[0]
    s