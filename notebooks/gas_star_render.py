import numpy as np
import pynbody
from pynbody.plot import sph
import matplotlib.pylab as plt
from util import np_printoptions
plt.style.use('dark_background')


def plot_gas_and_stars(sim, velocity_proj=False, **kwargs):
    """Wrapper around pynbody.plot.stars.render using correct smoothing length"""
    sim.properties.pop('boxsize', None)
    sim.g['smooth'] /= 2
    sim.s['smooth'] /= 2
    pynbody.analysis.halo.center(sim)
    if velocity_proj:
        """x-axis is aligned with the overall mean velocity of
        the snaphot and the vertical axis is the z axis rotated by the elevation angle
        of the velocity"""
        velocity = sim['vel'].mean(axis=0)
        backup = sim['pos'].copy()
        backup_v = sim['vel'].copy()
        v_xy = np.linalg.norm(velocity[0:2])
        v_x, v_y, v_z = velocity
        alpha = np.sign(v_y) * np.arccos(v_x/v_xy) * 180.0/np.pi
        theta = np.arctan(v_z/v_xy) * 180.0/np.pi
        r1=sim.rotate_z(alpha)
        r2=sim.rotate_y(theta)

    fig, (ax_g, ax_s) = plt.subplots(nrows=1, ncols=2, figsize=(16,8))
    try:
        snap = int(sim.filename[-4:])
#         with np_printoptions(precision=2):
#             title = '$t={:5.2f}$ Gyr, snap={}\nv = {}'.format(sim.properties['time'].in_units("Gyr"), snap, velocity)
        im = sph.image(sim.g, qty="rho", units="g cm^-2", subplot=ax_g, #title=title,
                       ret_im=True, show_cbar=False, **kwargs)
        width=kwargs.get("width", 20)
        rgbim = pynbody.plot.stars.render(sim, axes=ax_s, width=width, clear=False, plot=False, ret_im=True)
        ax_s.imshow(rgbim[::-1, :], extent=(-width / 2, width / 2, -width / 2, width / 2))
        ax_s.set_xlabel('x [' + str(sim.s['x'].units) + ']')
        ax_s.set_ylabel('y [' + str(sim.s['y'].units) + ']')
        ax_g.set_xlabel('x [' + str(sim.s['x'].units) + ']')
        ax_g.set_ylabel('y [' + str(sim.s['y'].units) + ']')
        
        fig.tight_layout()
        fig.subplots_adjust(top=0.92, bottom=0.15)
        cbar_ax = fig.add_axes([0.2,  0.05, 0.6, 0.02])
        fig.colorbar(im, cax=cbar_ax, orientation='horizontal').set_label("rho [g cm^-2]")
        
        title = '$t={:5.2f}$ Gyr, snap={}'.format(sim.properties['time'].in_units("Gyr"), snap)
        fig.suptitle(title)
    except Exception as e:
        raise(e)
    finally:
        sim.g['smooth'] *= 2
        sim.s['smooth'] *= 2
        if velocity_proj:
            # revert is costly (7% for each transformation w.r.t. the sph.image function)
            # and does not work when the transformation has been applied on a particle family
            # r2.revert()
            # r1.revert()
            sim['pos'] = backup
            sim['vel'] = backup_v
    return