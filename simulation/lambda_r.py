import warnings

import numpy as np


def lambda_r_formula(flux, r, v_los, sig_los):
    return np.sum(flux * r * np.abs(v_los)) / np.sum(flux * r * np.sqrt(v_los**2 + sig_los**2))


def lambda_r(flux, v_los, sig_los):
    resolution = flux.shape[-1]
    x = y = np.arange(-resolution/2, resolution/2)
    xx, yy = np.meshgrid(x, y)
    r = np.sqrt(xx**2 + yy**2)
    return np.sum(flux * r * np.abs(v_los)) / np.sum(flux * r * np.sqrt(v_los**2 + sig_los**2))


def lambda_r_r_lim(flux, v_los, sig_los, width, r_lim):
    """Compute lambda up to a certain radius `r_lim`"""
    assert flux.ndim == v_los.ndim == sig_los.ndim == 2
    resolution = flux.shape[-1]
    # print(f"flux.shape={flux.shape}")
    x = y = np.arange(-resolution/2, resolution/2)
    xx, yy = np.meshgrid(x, y)
    r = np.sqrt(xx**2 + yy**2)
    r_max = r_lim/width * resolution
    r[r>r_max] = 0
    return np.sum(flux * r * np.abs(v_los)) / np.sum(flux * r * np.sqrt(v_los**2 + sig_los**2))


def lambda_r_sim(flux, v_los, sig_los):
    assert flux.ndim == v_los.ndim == sig_los.ndim == 3
    resolution = flux.shape[-1]
    # print(f"flux.shape={flux.shape}")
    x = y = np.arange(-resolution/2, resolution/2)
    xx, yy = np.meshgrid(x, y)
    z = np.sqrt(xx**2 + yy**2)
    r = z[np.newaxis]
    # print(f"r.shape={r.shape}")
    return np.sum(flux * r * np.abs(v_los), axis=(1,2)) / np.sum(flux * r * np.sqrt(v_los**2 + sig_los**2), axis=(1,2))

def lambda_r_profile(flux, v_los, sig_los, width, r_lim, bins=50):
    resolution = flux.shape[-1]
    if bins>resolution/2:
        warnings.warn(f'bins > resolution/2 ({bins} > {resolution//2}): this does not make a lot of sense...')
    prof = list()
    r_arr = np.linspace(0, r_lim, bins)
    for r in r_arr:
        prof.append(lambda_r_r_lim(flux, v_los, sig_los, width=width, r_lim=r))
    return np.array(prof)


if __name__ == '__main__':
    # example usage
    import tqdm
    import matplotlib
    import matplotlib.pyplot as plt
    from simulation.simdata import get_maps
    sim_name = 'mb.71002_p100_a800_r600'
    maps = get_maps(sim_name, orbit_sideon=True)
    flux, v_los, sig_los = maps['lum'], maps['vlos'], maps['sig_los']
    n_prof = len(flux)
    bins = 50
    width = maps.meta['WIDTH']
    norm = matplotlib.colors.Normalize(vmin=0, vmax=n_prof)
    cmap = matplotlib.cm.get_cmap('jet', 12)
    mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    fig, ax = plt.subplots()
    x = np.linspace(0, width/2, bins)
    prof_list = list()
    for i in tqdm.tqdm(range(len(flux))):
        prof_list.append(lambda_r_profile(flux[i], v_los[i], sig_los[i], width))

    for i in range(0, len(flux), 10):
        ax.plot(x, prof_list[i], color=mappable.to_rgba(i))
    ax.set_ylabel('$\lambda_R$')
    ax.set_xlabel('r [kpc]');
    fig.colorbar(mappable)
    plt.show()