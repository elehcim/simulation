import os
import pynbody
from pynbody.plot import sph
import matplotlib.pylab as plt

from notebooks.snap_io import load_moria_sim_and_kicked, load_moria, load_kicked, load_sim

SIMNUMBER = "69002_p200.0_a600.0_r600.0_c8.15"
# SIMNUMBER = "69002_p598.0_a600.0_r600.0_c8.15"
kicked = True

faceon = False

resolution = 500
width = 20 # kpc
vmax=2e-1
vmin=5e-4
figsize = (15,15)

plt.style.use("dark_background")


import contextlib
@contextlib.contextmanager
def np_printoptions(*args, **kwargs):
    original = np.get_printoptions()
    np.set_printoptions(*args, **kwargs)
    try:
        yield
    finally: 
        np.set_printoptions(**original)

def gas_image(sim, **kwargs):
    sim.properties.pop('boxsize', None)
    sim.g['smooth'] /= 2
    pynbody.analysis.halo.center(sim, shrink_factor=0.9)
    if faceon:
        pynbody.analysis.angmom.faceon(sim)
    try:
        img = sph.image(sim.g, qty="rho", units="g cm^-2", **kwargs)
    finally:
        sim.g['smooth'] *= 2
    return img

def gas_velocity_lateral_view(sim, **kwargs):
    sim.properties.pop('boxsize', None)
    sim.g['smooth'] /= 2
    pynbody.analysis.halo.center(sim, shrink_factor=0.9)
    velocity = sim['vel'].mean(axis=0)
    backup = sim['pos'].copy()
    v_modulus = np.linalg.norm(velocity)
    x_p = (velocity[1]  * sim['pos'][:,0] + velocity[0] * sim['pos'][:,1])/v_modulus 
    y_p = (-velocity[0] * sim['pos'][:,0] + velocity[1] * sim['pos'][:,1])/v_modulus
    sim['pos'][:,0] = x_p
    sim['pos'][:,1] = y_p
    try:
        snap = int(sim.filename[-4:])
        with np_printoptions(precision=2):
            title = '$t={:5.2f}$ Gyr, snap={}\nv = {} kpc/s'.format(sim.properties['time'].in_units("Gyr"), snap, velocity)
        img = sph.image(sim.g, qty="rho", units="g cm^-2", title=title, **kwargs)
    finally:
        sim.g['smooth'] *= 2
        sim['pos'] = backup
    return img

def gas_orthogonal_projection(sim, qty='rho', **kwargs):
    sim.properties.pop('boxsize', None)
    sim.g['smooth'] /= 2
    pynbody.analysis.halo.center(sim, shrink_factor=0.9)
    fig, ((ax_xy, ax_zy), (ax_xz, ax_angmom)) = plt.subplots(nrows=2, ncols=2, figsize=(15,15))
    kwargs['ret_im']=True
    if "units" not in kwargs:
        kwargs['units'] = "g cm^-2"
    u_st = sim['pos'].units.latex()
    backup = sim['pos'].copy()
    try:
        im = sph.image(sim.g, subplot=ax_xy, **kwargs)
        ax_xy.set_xlabel("$x/%s$" % u_st)
        ax_xy.set_ylabel("$y/%s$" % u_st)
        
        x=sim['pos'][:,0].copy()
        sim['pos'][:,0] = sim['pos'][:,2].copy()
        sph.image(sim.g, subplot=ax_zy, **kwargs)
        ax_zy.set_xlabel("$z/%s$" % u_st)
        ax_zy.set_ylabel("$y/%s$" % u_st)

        sim['pos'][:,1] = sim['pos'][:,2]
        sim['pos'][:,0] = x
        sph.image(sim.g, subplot=ax_xz, **kwargs)
        ax_xz.set_xlabel("$x/%s$" % u_st)
        ax_xz.set_ylabel("$z/%s$" % u_st)

        pynbody.analysis.angmom.faceon(sim)
        sim['pos'] = backup
        sph.image(sim.g, subplot=ax_angmom, **kwargs)
        ax_angmom.set_xlabel("$x/%s$" % u_st)
        ax_angmom.set_ylabel("$y/%s$" % u_st)
        ax_angmom.set_title("angular momentum faceon")
        
        fig.tight_layout()
        fig.subplots_adjust(top=0.96, bottom=0.11)
        cbar_ax = fig.add_axes([0.2,  0.05, 0.6, 0.02])
        fig.colorbar(im, cax=cbar_ax, orientation='horizontal').set_label(qty+"/"+kwargs['units'])
        
        snap = int(sim.filename[-4:])
        title = '$t={:5.2f}$ Gyr, snap={}'.format(sim.properties['time'].in_units("Gyr"), snap)
        fig.suptitle(title)
    finally:
        sim.g['smooth'] *= 2
        sim['pos'] = backup
    return fig


def main():
    snap_list = load_kicked(SIMNUMBER) if kicked else load_moria(SIMNUMBER)
    
    folder = "pngs_{}".format(SIMNUMBER)
    os.makedirs(folder, exist_ok=True )
    
    for sim in snap_list:
        snap = int(sim.filename[-4:])
        fig, ax = plt.subplots(figsize=figsize)
        if not faceon:
            filename = os.path.join(folder,"gas_image_{:03d}.png".format(snap))
            title = '$t={:5.2f}$ Gyr, snap={:03d}'.format(sim.properties['time'].in_units("Gyr"), snap)
        else:
            filename = os.path.join(folder,"gas_image_angmomfaceon_{:03d}.png".format(snap))
            title = '$t={:5.2f}$ Gyr, snap={:03d} (angmomfaceon)'.format(sim.properties['time'].in_units("Gyr"), snap)
        gas_image(sim, width=width, resolution=resolution, vmin=vmin, vmax=vmax, title=title, filename=filename);
        print("Saved", filename)

def save_orthogonal():
    snap_list = load_kicked(SIMNUMBER) if kicked else load_moria(SIMNUMBER)
    
    folder = "pngs_{}".format(SIMNUMBER)
    os.makedirs(folder, exist_ok=True )
    
    for sim in snap_list:
        snap = int(sim.filename[-4:])
        filename = os.path.join(folder,"gas_image_proj_{:03d}.png".format(snap))
        fig = gas_orthogonal_projection(sim, width=width, resolution=resolution, vmin=vmin, vmax=vmax);
        fig.savefig(filename)
        print("Saved", filename)

def save_velocity_lateral_view():
    snap_list = load_kicked(SIMNUMBER) if kicked else load_moria(SIMNUMBER)
    
    folder = "pngs_{}".format(SIMNUMBER)
    os.makedirs(folder, exist_ok=True )
    
    for sim in snap_list:
        snap = int(sim.filename[-4:])
        fig, ax = plt.subplots(figsize=figsize)
        filename = os.path.join(folder,"gas_image_vel_proj_{:03d}.png".format(snap))
        title = '$t={:5.2f}$ Gyr, snap={:03d}'.format(sim.properties['time'].in_units("Gyr"), snap)
        gas_image(sim, width=width, resolution=resolution, vmin=vmin, vmax=vmax, title=title, filename=filename);
        print("Saved", filename)

if __name__ == '__main__':
    # main()
    # save_orthogonal()
    save_velocity_lateral_view()
