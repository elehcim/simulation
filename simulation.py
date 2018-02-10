import numpy as np
import matplotlib.pylab as plt
import logging
import pynbody
from snap_io import load_moria_sim_and_kicked, load_moria, load_kicked, load_sim
from util import np_printoptions

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)

def mass_resolution(snap):
    return (snap['mass'].sum()/len(snap)).in_units("Msol")

def velocity_projection(sim):
    v_x, v_y, v_z = sim['vel'].mean(axis=0)
    v_xy = np.linalg.norm(velocity[0:2])
    alpha = np.sign(v_y) * np.arccos(v_x/v_xy) * 180.0/np.pi
    theta = np.arctan(v_z/v_xy) * 180.0/np.pi            
    return alpha, theta


class Simulation(object):
    """docstring for Simulation"""
    times = None
    cog = None
    def __init__(self, sim_id):
        self.sim_id = sim_id
        self.snap_list = load_sim(sim_id)

    def snap(self, idx):
        return self.snap_list[idx]

    def get_times(self):
        self.times = np.zeros(len(self), dtype=float)
        for i, snap in enumerate(self.snap_list):
            self.times[i] = snap.properties['time'].in_units('Gyr')

    def mass_resolution(self):
        pass

    def __getitem__(self, idx):
        return self.snap_list(idx)

    def __len__(self):
        return len(self.snap_list)

    def __repr__(self):
        return "{}: ({}) {}".format(self.sim_id, len(self), self.snap(0).__repr__())
    
    def compute_cog(self, save_cache=False, cache_file=None, verbose=True, family=None):
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
        self.cog = np.zeros((3, len(self)), dtype=float)
        for i, snap in enumerate(snapshots):
            if family is not None:
                snap = snap.__getattr__(family)

            if verbose:
                print("{:03d} Analysing {} (time {:.4f} Gyr)".format(i, snap.filename, snap.properties['time'].in_units('Gyr')))
            
            mass = snap['mass']
            pos = snap['pos']
            tot_mass = mass.sum()
            self.cog[:,i] = np.sum(mass * pos.transpose(), axis=1) / tot_mass
            
            if save_cache:
                if cache_file is not None:
                    cache_file = self.sim_id + ".cog.npz"
                np.savez(cache_file, cog=cog)
        self.cog = cog
        return cog

    def _center_all(self):
        for snap in self.snap_list:
            pynbody.analysis.halo.center(snap)

class MoriaSim(Simulation):
    """docstring for MoriaSim"""
    def __init__(self, sim_id, kicked=False):
        self.sim_id = sim_id
        logger.info("loading simulation: {}".format(sim_id))
        self.snap_list = load_kicked(sim_id) if kicked else load_moria(sim_id)
        self.kicked = kicked
        # Remove boxsize which complicates the plotting
        for i, snap in enumerate(self.snap_list):
            if i==0:
                self.boxsize = snap.properties.pop('boxsize', None).copy()
            snap.properties.pop('boxsize', None)
    
    def clear_cache(self):
        pass

    @property
    def t_range(self):
        if self.times is None:
            self.get_times()
        return self.times.min(), self.times.max()

    def plot_star(self, snap):
        """Wrapper around pynbody.plot.stars.render using correct smoothing length"""
        pass    

    def plot_gas(self, i, faceon=False, **kwargs):
        snap = self.snap_list[i]
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


    def plot_gas_and_stars(self, i, velocity_proj=False, sfh=False, **kwargs):
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
            alpha, theta = velocity_projection(snap)
            r1=snap.rotate_z(alpha)
            r2=snap.rotate_y(theta)

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
            rgbim = pynbody.plot.stars.render(snap, axes=ax_s, width=width, clear=False, plot=False, ret_im=True)
            ax_s.imshow(rgbim[::-1, :], extent=(-width / 2, width / 2, -width / 2, width / 2))
            ax_s.set_xlabel('x [' + str(snap.s['x'].units) + ']')
            ax_s.set_ylabel('y [' + str(snap.s['y'].units) + ']')
            ax_g.set_xlabel('x [' + str(snap.s['x'].units) + ']')
            ax_g.set_ylabel('y [' + str(snap.s['y'].units) + ']')
            
            fig.tight_layout() # only plots above are affected
            fig.subplots_adjust(top=0.92, bottom=0.15)
            cbar_ax = fig.add_axes([0.2,  0.06, 0.6, 0.02])
            fig.colorbar(im, cax=cbar_ax, orientation='horizontal').set_label("rho [g cm^-2]")
            if sfh:
                #  [left, bottom, width, height]
                ax_sfh = fig.add_axes([0.1,  -0.3, 0.35, 0.26])
                # ignore AccuracyWarning that is issued when an integral is zero
                import warnings
                from scipy.integrate.quadrature import AccuracyWarning
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", category=AccuracyWarning)
                    my_range = (0, 13.7)
                    pynbody.plot.stars.sfh(snap, trange=my_range, range=my_range, subplot=ax_sfh)
                ax_sfh.axvline(x=snap_time_gyr, linestyle="--")
                ax_sfh.set_xlabel("Time [Gyr]")
                ax_sfh.set_ylabel("SFR [M$_\odot$ yr$^{-1}$]")
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
        return


    def interact_plot_gas_star():
        pass





def time_range_kicked_moria():
    trange = (min(np.min(times_moria), np.min(times_kicked)), max(np.max(times_moria), np.max(times_kicked)))
    trange = (0, min(np.max(times_moria), np.max(times_kicked)))
    bins = np.linspace(*trange, bins_sfr)
    trange, times_moria[-1], times_kicked[-1]

if __name__ == '__main__':
    SIMNUMBER = "69002_p200.0_a600.0_r600.0_c8.15"
    kicked=True
    snap_list = load_kicked(SIMNUMBER) if kicked else load_moria(SIMNUMBER)
    sim = MoriaSim(SIMNUMBER, kicked)
    s = snap_list[0]
    s