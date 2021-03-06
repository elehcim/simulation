import csv
import logging
import pickle
import os
from functools import lru_cache
from multiprocessing import Process, Queue

import matplotlib.pylab as plt
import numpy as np
import pandas as pd
import pynbody
from .analyze_sumfiles import get_sumfile
from .parsers.parse_trace import parse_trace, parse_dens_trace
from .parsers.parse_info import parse_info
from .snap_io import load_moria, load_kicked, load_sim, make_snaps_path, snapshot_file_list
from .util import np_printoptions, make_df_monotonic_again_using_reference_df, get_sim_name
from .plot.plot_trace import plot_trace, plot_trace_df
from .units import gadget_time_units, gadget_dens_units, gadget_vel_units
from .simdata import SIM_NAME_DICT, SIMS_DIR

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
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


def _my_cog(snap):
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


def mu(temp):
    """Mean atomic mass"""
    x = np.empty(len(temp)).view(pynbody.array.SimArray)
    x[np.where(temp >= 1e4)[0]] = 0.59
    x[np.where(temp < 1e4)[0]] = 1.3
    x.units = pynbody.units.Unit("1")
    return x


def cs(temp):
    """Sound speed"""
    _mu = mu(temp)
    return np.sqrt(5.0 / 3.0 * pynbody.units.k * temp / _mu / pynbody.units.m_p).in_units('km s**-1')


def get_param_used(path):
    """
    Parse the `parameters-usedvalues` files

    Parameters
    ----------
    path : str
        path of the simulation `out` directory

    Returns
    -------
    A dictionary of {parameter: value}. Values are all strings.
    """
    d = {}
    path = os.path.expanduser(path)
    if os.path.isdir(path):
        path = os.path.join(path, 'parameters-usedvalues')
    else:
        return None
    try:
        with open(path) as csvfile:
            logger.debug("Found parameter file")
            spamreader = csv.reader(csvfile, delimiter=' ', skipinitialspace=True)
            for line in spamreader:
                try:
                    k, v = line
                    d[k] = v
                except ValueError as e:
                    logger.warning("{}. But continuing".format(e))
    except FileNotFoundError as e:
        logger.warning("Parameter file not found: {}. But continuing".format(e))
        return None
    return d


def get_dens_trace(path):
    path = os.path.expanduser(path)
    if os.path.isdir(path):
        path = os.path.join(path, 'dens_temp_trace.txt')
    else:
        return None
    try:
        df = parse_dens_trace(path)
        logger.info("Found dens_temp_trace file")
    except FileNotFoundError as e:
        logger.warning("dens_temp_trace file not found: {}. But continuing".format(e))
        return None
    return df


def get_trace(path):
    path = os.path.expanduser(path)
    if os.path.isdir(path):
        if os.path.isfile(os.path.join(path, 'trace.pkl')):
            logger.info("Found cached trace file")
            return pd.read_pickle(os.path.join(path, 'trace.pkl'))
        else:
            path = os.path.join(path, 'trace.txt')
    else:
        return None

    try:
        df = parse_trace(path)
        logger.info("Found trace file")
    except FileNotFoundError as e:
        logger.warning("trace file not found: {}. But continuing".format(e))
        return None
    return df

def get_info(path):
    path = os.path.expanduser(path)
    if os.path.isdir(path):
        if os.path.isfile(os.path.join(path, 'info.pkl')):
            logger.info("Found cached info file")
            return pd.read_pickle(os.path.join(path, 'info.pkl'))
        else:
            path = os.path.join(path, 'info.txt')
    else:
        return None

    try:
        df = parse_info(path)
        logger.info("Found info file")
    except FileNotFoundError as e:
        logger.warning("info file not found: {}. But continuing".format(e))
        return None
    return df

def get_compiler_options(path):
    path = os.path.expanduser(path)
    if os.path.isdir(path):
        path = os.path.join(path, 'compiler.txt')
    else:
        return None
    try:
        with open(path) as f:
            logger.info("Found compiler file")
            ll = [s.strip() for s in f.readlines()]
    except FileNotFoundError:
        return None
    return ll

def get_times_header(sim_dir):
    """Faster way to get times in the header"""
    snap_name_list = snapshot_file_list(os.path.expanduser(sim_dir), include_dir=True)
    return np.array([pynbody.load(snap).header.time for snap in snap_name_list])


class Simulation:
    """An object to work with Gadget2 simulations.

    The usual content of a simulation folder which `Simulation` can read is the following:

    simdir:
      |- Gadget2 (executable)
      |- Makefile
      |- out/
        |- compiler.txt
        |- cpu.txt
        |- energy.txt
        |- info.txt
        |- parameter-usedvalues
        |- sfc.dat
        |- timings.txt
        |- [trace.txt]
        |- [dens_temp_trace.txt]
    """
    _times = None
    cog = None
    _computed_cog = False

    def __init__(self, sim_dir, snap_indexes=None, sim_id=None, force_cosmo=False):
        """
        Initialize a Simulation.

        Parameters
        ----------
        sim_dir : str
            path of the simulation `out` dir
        sim_id : str
            Identification ID of the simulation. If None, `sim_dir` is used
        snap_indexes : iterable (list, slice, tuple, int)
            list of the snapshots to load. Example to take the last 10 snaps:
                        `snap_indexes = slice(-10,None)`
        force_cosmo : bool
            Use default cosmology, do not read cosmological parameters from the snapshots.
            It can be useful if the snapshots contain an incorrect value of OmegaM0 for example.
        """
        # TODO join this __init__ with the Moria or Kicked class

        # Shortcut
        if sim_dir in SIM_NAME_DICT:
            short_name = sim_dir
            sim_dir = os.path.join(SIMS_DIR, SIM_NAME_DICT[sim_dir], 'out')
        else:
            short_name = None

        if sim_id is None:
            sim_id = sim_dir
        self.sim_id = sim_id
        self._sim_dir = sim_dir
        logger.info("loading simulation: {}".format(sim_id))
        self.params = get_param_used(sim_dir)
        self.compiler_opts = get_compiler_options(sim_dir)
        self._snap_indexes = snap_indexes if snap_indexes is not None else slice(None)
        self.snap_list = self._load(sim_dir, force_cosmo, snap_indexes)
        if len(self.snap_list) == 0:
            raise RuntimeError("No snaphots found in {}".format(sim_dir))

        self._centered = np.zeros(len(self.snap_list), dtype=bool)

        # I do this first because I need it (non monotonic version) for the trace.txt
        self.dens_trace = get_dens_trace(sim_dir)
        if self.dens_trace is not None:
            # np.digitize works only if it is monotonic
            if not self.dens_trace.t.is_monotonic_increasing:
                logger.info("dens_temp_trace.txt file is non-monotonic. Trying to recover")
                df_mono = make_df_monotonic_again_using_reference_df(self.dens_trace, self.dens_trace)
                self._dens_trace_orig = self.dens_trace.copy()
                self.dens_trace = df_mono

        self.trace = get_trace(sim_dir)
        if self.trace is not None:
            # np.digitize works only if t is monotonic
            if not self.trace.t.is_monotonic_increasing:
                logger.info("trace.txt file is non-monotonic. Trying to recover")
                # I assume that if this is non-monotonic also dens_trace is, so that I have a _dens_trace_orig version
                df_mono = make_df_monotonic_again_using_reference_df(self.trace, self._dens_trace_orig)
                self._trace_orig = self.trace.copy()
                self.trace = df_mono
            assert self.trace.t.is_monotonic_increasing, "trace is still non monotonically increasing"

    def _load(self, sim_id, force_cosmo=False, snap_indexes=None):
        snap_name_list = snapshot_file_list(os.path.expanduser(sim_id), include_dir=True)
        logger.info("Found {} snapshots".format(len(snap_name_list)))
        if snap_indexes is not None:
            if isinstance(snap_indexes, (list, tuple)):
                snap_name_list = np.array(snap_name_list)[snap_indexes]
            else:
                snap_name_list = snap_name_list[snap_indexes]
            logger.info("Taking {} snapshots ({})".format(len(snap_name_list), snap_indexes))
        if pynbody.__version__ >= '0.47':
            snap_list = list(pynbody.load(snap, ignore_cosmo=True) for snap in snap_name_list)
        else:
            snap_list = list(pynbody.load(snap) for snap in snap_name_list)

        logger.info('Loading cosmological parameters')
        for i, snap in enumerate(snap_list):
            if self.params is not None:
                snap.properties['h'] = float(self.params['HubbleParam'])
                snap.properties['omegaL0'] = float(self.params['OmegaLambda'])
                snap.properties['omegaM0'] = float(self.params['OmegaBaryon'])
            # Take boxsize from first snap
            if i == 0:
                self.boxsize = snap.properties.get('boxsize', None)
            if force_cosmo:
                if i == 0:
                    logger.info('Forcing cosmological parameters (h=0.7, omegaL0=0.72, omegaM0=0.28)')
                snap.properties['h']= 0.7
                snap.properties['omegaL0']= 0.72
                snap.properties['omegaM0']= 0.28
        return snap_list

    def _cog(self, i):
        snap = self.snap_list[i]
        mass = snap['mass']
        pos = snap['pos']
        tot_mass = mass.sum()
        return np.sum(mass * pos.transpose(), axis=1) / tot_mass

    def _save_trace_cache(self):
        trace_name = 'trace.pkl'
        logger.info("Saving  cached trace file")
        self.trace.to_pickle(os.path.join(self._sim_dir, trace_name))

    def snap(self, idx):
        return self.snap_list[idx]

    def get_traj(self):
        if not self.trace.t.is_monotonic_increasing:
            logger.error("trace.txt file is non-monotonic. Cannot get trajectory")
            raise RuntimeError("trace.txt file is non-monotonic. Cannot get trajectory")
        # self.trace should be monotonic here
        locations = np.digitize(self.times.in_units(gadget_time_units), self.trace.t, right=True)
        x = self.trace.x[locations]
        y = self.trace.y[locations]
        z = self.trace.z[locations]
        return pynbody.array.SimArray(np.vstack([x, y, z]), 'kpc').T

    def get_rv_host(self):
        if self.dens_trace is not None:
            locations = np.digitize(self.times.in_units(gadget_time_units), self.dens_trace.t, right=True)
            rho_host  = pynbody.array.SimArray(self.dens_trace.rho[locations], gadget_dens_units)
            v_host = pynbody.array.SimArray(self.dens_trace.vel[locations],  gadget_vel_units)
        else:
            rho_host = v_host = None
        return rho_host, v_host

    def get_temp_host(self):
        if self.dens_trace is not None:
            locations = np.digitize(self.times.in_units(gadget_time_units), self.dens_trace.t, right=True)
            temp = pynbody.array.SimArray(self.dens_trace.temp[locations], 'K')
        else:
            temp = None
        return temp

    def cs(self):
        temp = self.get_temp_host()
        return cs(temp)

    def mach(self):
        _, v = self.get_rv_host()
        return v / self.cs()

    # def get_omega(self):
    #     omega = None
    #     if self.is_moving_box:
    #         omega_arr = get_omega_box(self)
    #         omega = pynbody.array.SimArray(omega_arr, 1/gadget_time_units)
    #     return omega

    @property
    def peri(self):
        if self.is_moving_box:
            sim_name = get_sim_name(self._sim_dir)
            i = sim_name.find('_p')
            j = sim_name.find('_a')
            return sim_name[i+2:j]
        else:
            return ""

    @property
    def is_moving_box(self):
        return self.trace is not None
        # This should be the correct way but it is not usable at the moment:
        # given that up to now the simulations do not save the HOT_HALO flag on file
        # if self.compiler_opts is not None:
        #     return 'HOT_HALO' in self.compiler_opts
        # else:
        #     return False

    @property
    def properties(self):
        d = self.snap_list[0].properties.copy()
        for k in ('time', 'a'):
            if k in d:
                del d[k]
        # del d['time'], d['a']
        d['time_begin'] = self.snap_list[0].properties['time'].in_units('Gyr')
        d['time_end'] = self.snap_list[-1].properties['time'].in_units('Gyr')
        d['z_begin'] = self.snap_list[0].header.redshift
        d['z_end'] = self.snap_list[-1].header.redshift
        return d

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
        return pynbody.array.SimArray([snap.properties['time'].in_units('Gyr') for snap in self.snap_list], units=pynbody.units.Gyr)

    @property
    @lru_cache(1)
    def times_header(self):
        return np.array([snap.header.time for snap in self.snap_list])

    @property
    @lru_cache(1)
    def ram_pressure(self):
        return self.dens_trace.vel**2 * self.dens_trace.rho

    @property
    def r(self):
        if self.is_moving_box:
            locations = np.digitize(self.times.in_units(gadget_time_units), self.trace.t, right=True)
            # self.r = pynbody.array.SimArray(self.trace.r[locations], 'kpc')
            return self.trace.r[locations]
        elif self.cog is not None:
            return np.linalg.norm(self.cog, axis=0)

    def get_times(self):
        self._times = np.zeros(len(self))
        for i, snap in enumerate(self.snap_list):
            self._times[i] = snap.properties['time'].in_units('Gyr')
        return self._times

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

    # def is_id_duplicated(self):
    #     # FIXME
    #     """For each snapshot return True if some IDs (iord) are duplicated"""
    #     is_dup = [contains_duplicated(snap['iord']) for snap in self.snap_list]
    #     return np.array(is_dup, dtype=bool)

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
            save_cache: Save results in a npz file
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
            cache_file = os.path.join(self._sim_dir, cache_dir,
                os.path.basename(self.sim_id) + (".cog.npz" if family is None else ".{}.cog.npz".format(family)))
            if save_cache:
                os.makedirs(os.path.dirname(cache_file), exist_ok=True)

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

    #TODO
    """
    def m_star(self, ax_m_star, snap_time_gyr=None, last_snap=-1, **kwargs):
        if ax_m_star is None:
            fig, ax_m_star = plt.subplots(1, figsize=(8,6))
        # ignore AccuracyWarning that is issued when an integral is zero
        import warnings
        from scipy.integrate.quadrature import AccuracyWarning
        np.zeros(len(self.snap_list), dtype=float)
        self.snap_list[last_snap]
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=AccuracyWarning)
            pynbody.plot.stars.m_star(self.snap_list[last_snap], subplot=ax_m_star, **kwargs)
        if snap_time_gyr is not None:
            ax_m_star.axvline(x=snap_time_gyr, linestyle="--")
        # ax_m_star.set_title("m_star")
        ax_m_star.set_xlabel("Time [Gyr]")
        ax_m_star.set_ylabel("SFR [M$_\odot$ yr$^{-1}$]")
        return ax_m_star
    """

    def plot_gas_and_stars(self, i, velocity_proj=False, sfh=False, cog=False, starsize=None, **kwargs):
        from .sfh_in_box import plot_binned_sfh
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
                if self.is_moving_box:
                    plot_binned_sfh(self, bins=100, ax=ax_sfh, drawstyle='steps')
                    ax_sfh.axvline(x=snap_time_gyr, linestyle="--")
                else:
                    self.plot_sfh(ax_sfh, snap_time_gyr)

            if cog:
                ax_cog = fig.add_axes([0.63,  -0.3, 0.26, 0.26])
                if self.is_moving_box:
                    plot_trace_df(self.trace, ax=ax_cog)
                    idx = int(np.digitize(snap.properties['time'], self.trace.t))
                    cur_pos = self.trace.iloc[idx]
                    ax_cog.scatter(cur_pos.x, cur_pos.y, color="red", marker='o', facecolors='none')
                    ax_cog.scatter(0, 0, marker='+', color="b")
                else:
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

    def __init__(self, sim_id, kicked=False, snap_indexes=None):
        self.sim_id = str(sim_id)
        self.kicked = kicked

        # self.snap_list = load_kicked(sim_id) if kicked else load_moria(sim_id)
        # super(MoriaSim, self).__init__(sim_id)
        self._load(snap_indexes=snap_indexes)
        sumfile_path = os.path.join(self._sf_moria, self.sim_id + ".dat")
        if os.path.isfile(sumfile_path):
            logger.info("Getting sumfile: {}".format(sumfile_path))
            self.sumfile = get_sumfile(os.path.join(self._sf_moria, self.sim_id + ".dat"))
        else:
            logger.info("No sumfile found")
        self._centered = np.zeros(len(self.snap_list), dtype=bool)

    def _load(self, snap_indexes=None):
        logger.info("loading simulation: {}".format(self.sim_id))
        _loader = load_kicked if self.kicked else load_moria
        self.snap_list = _loader(self.sim_id, snap_indexes=snap_indexes)
        self._sim_dir = make_snaps_path(self.sim_id, self.kicked)

        # Overwrite and fix cosmological parameters
        print('Fixing cosmological parameters of MoRIA simulation')
        for i, snap in enumerate(self.snap_list):
            snap.properties['h']= 0.7
            snap.properties['omegaL0']= 0.72
            snap.properties['omegaM0']= 0.28
            # Remove boxsize which complicates the plotting
            if i==0 and 'boxsize' in snap.properties:
                self.boxsize = snap.properties.pop('boxsize').copy()
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


# class BhSim(Simulation):
#     """docstring for MoriaSim"""
#     _sim_dir = "/home/michele/sim/MySimulations/bh"
#     _sf_bh_folder = os.path.join(_sim_dir, "results/sumfiles/")
#
#     def __init__(self, sim_id):
#         self.sim_id = sim_id
#         # self.snap_list = load_kicked(sim_id) if kicked else load_moria(sim_id)
#         # super(MoriaSim, self).__init__(sim_id)
#         self._load(sim_id)
#
#     def _load(self, sim_id):
#         logger.info("loading simulation: {}".format(sim_id))
#         self.snap_list = load_sim(os.path.join(self._sim_dir, sim_id, "out"))
#         # Remove boxsize which complicates the plotting
#         for i, snap in enumerate(self.snap_list):
#             if i==0:
#                 self.boxsize = snap.properties.pop('boxsize', None).copy()
#             snap.properties.pop('boxsize', None)

if __name__ == '__main__':
    SIMNUMBER = "69002_p200.0_a600.0_r600.0_c8.15"
    kicked=True
    sim = MoriaSim(SIMNUMBER, kicked)
    s = sim.snap_list[0]
