fig, axes = plt.subplots(ncols=4, nrows=2, figsize=(22, 10), dpi=300)
[[ax_lam, ax_ellip, ax_reff, ax_sfh],[ax_mu, ax_sigma, ax_n, ax_orb]] = axes
values_to_plot = ('lambda_r', 'ellipticity', 'r_eff', 'sfr', 'mu_e', 'sigma', 'n')

axs = dict(zip(values_to_plot, axes.flatten()[:-1]))

axs['lambda_r'].plot(df.t, df.lambda_r_mean)
# axs['lambda_r'].scatter(df.time, df.lambda_r, color='g', alpha=0.4)
axs['lambda_r'].set_ylabel('$\lambda_R$')
axs['lambda_r'].set_title('Specific Stellar Angular Momentum')
# axs['lambda_r'].set_ylim(0, 1)

axs['ellipticity'].plot(df.t, df.ellipticity_mean)
axs['ellipticity'].set_ylabel('$\epsilon$')
axs['ellipticity'].set_title("Ellipticity")
# axs['ellipticity'].set_ylim(0, 0.8)

# axs['reff'].plot(df.time, df.r_eff)
axs['reff'].plot(df.t, df.r_eff3d_mean)
# axs['reff'].plot(df.t, df.r_eff_mean)
# axs['reff'].plot(df.time, df.r_eff_fit)
axs['reff'].set_ylabel('$R_e$ [kpc]')
axs['reff'].set_title("Effective radius")
# axs['reff'].set_ylim(0, None)

axs['sfr'].plot(df.t, df.sfr_mean)
axs['sfr'].set_ylabel('SFR [Msol/yr]')
axs['sfr'].set_title("SFH")

# axs['tot_mag'].plot(df.t, df.mag_v_mean)
# axs['tot_mag'].set_ylabel('$m_V$')
# axs['tot_mag'].invert_yaxis()
# axs['tot_mag'].set_title('Total Magnitude (V band)')

# axs['tot_mag'].plot(df.t, df.metals_star)
# axs['tot_mag'].set_ylabel('$[Fe/H]_\star$')
# axs['tot_mag'].set_title('Metallicity of stars')

axs['sigma'].plot(df.t, df.sigma_star_mean, label='$\sigma_\star$')
# axs['sigma'].plot(df.t, df.sigma_gas_mean, label='$\sigma_g$')
axs['sigma'].set_ylabel('$\sigma [km/s]$')
axs['sigma'].set_title('$\sigma$')
axs['sigma'].legend()
# axs['sigma'].set_ylim(0, None)

axs['mu_e'].plot(df.t, df.mu_e_mean)
axs['mu_e'].invert_yaxis()
axs['mu_e'].set_ylabel('$\mu_e$ [mag/arcsec$^2$]')
axs['mu_e'].set_title('Surface Brightness at $R_e$')

axs['n'].plot(df.t, df.n_mean)
# axs['n'].plot(df.time, df.n_mean)
axs['n'].set_ylabel('n')
# axs['n'].set_ylim(0, None)
axs['n'].set_title('Sersic Index')

ax_orb.plot(df.x, df.y)
ax_orb.plot((0,0), 'r+')
ax_orb.set_aspect('equal')
ax_orb.set_xlim(-900, 900)
ax_orb.set_ylim(-900, 900)
ax_orb.set_xlabel('x [kpc]')
ax_orb.set_ylabel('y [kpc]')
ax_orb.set_title('Orbit')

plot_list = ('lambda_r_mean', 'ellipticity_mean', ('r_eff3d_mean'), 'sfr_mean', 'mu_e_mean', ('sigma_star_mean'), 'n_mean')

def fill_between(df, plottable):
    ser = df[plottable[:-5]]
    ser_mean = df[plottable]
    ser_std = df[plottable[:-5]+'_std']
    return ser_mean - ser_std, ser_mean + ser_std

for i, ax in enumerate(axes.flatten()[:-1]):
    ax.set_xlabel('Time [Gyr]')
    ax2 = ax.twinx()
    ax2.plot(df.t, df.r, "r--", alpha=0.2)
    ax2.set_yticks([200, 500, 800])
    ax2.tick_params(axis='both', direction='in')
    ax.tick_params(axis='both', direction='in')
    ax.set_xlabel('Time [Gyr]')
    ax.grid()

    plottable = plot_list[i]
    if isinstance(plottable, tuple) and len(plottable) > 1:
        for p in plottable:
            ax.plot(df_m.t, df_m[p])
    else:
        ax.plot(df_m.t, df_m[plottable])

    color_sim = ax.lines[0].get_color()
    color_moria = ax.lines[1].get_color()
    alpha = 0.4

#     ax.scatter(df.t, df[plottable[:-5]], s=1, color=color_sim, alpha=alpha)
#     ax.scatter(df_m.t, df_m[plottable[:-5]], s=1, color=color_moria, alpha=alpha)

    ax.fill_between(df.t, *fill_between(df, plottable), color=color_sim, alpha=0.2)
    ax.fill_between(df_m.t, *fill_between(df_m, plottable), color=color_moria, alpha=0.2)

#     sns.regplot(x="t", y=plottable[:-5], data=df, ax=ax, x_estimator=np.mean, truncate=True)

    ax.set_xlim(7.8, None)
    if i != 4:
        ax.set_ylim(0, None)

# ax_mu.set_ylim(None, None)

plt.subplots_adjust(wspace=0.3, hspace=0.2)
fig.suptitle(sim_name);
fig.savefig(sim_name + '_overview.png')