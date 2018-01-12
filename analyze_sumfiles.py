import astropy.units as u
from astropy.io import ascii
from astropy.table import Table

header = ["nr","time","redshift","ngas","nstar","ndark","Mgas","Mstar","Mdark","MHI","r_e_L","r_e_M","r_e_M_DM",
          "flatStar","flatDM","velDisp_star","velDisp_DM","kinE_star","kinE_DM","bound_star","bound_DM","SFR",
          "M_U","M_B","M_V","M_R","M_I","M_J","M_H","M_K","metall_star","FeH_star","MgFe_star","metall_lum","FeH_lum",
          "age_lum","Mgas_re","Mstar_re","Mdark_re","rcom_x","rcom_y","rcom_z","vcom_x","vcom_y","vcom_z",
          "L_x","L_y","L_z","I0","sig(I0)","r0","sig(r0)","n","sig(n)","mu0","sig(mu0)",
          "gas<1Re","gas<2Re","gas<3Re","gas<4Re","gas<5Re","gas<10Re","gas<30Re","gas_all",
          "DM<1Re","DM<2Re","DM<3Re","DM<4Re","DM<5Re","DM<10Re","DM<30Re","DM_all"]

_dl = u.dimensionless_unscaled
units = {"time":u.kpc/u.km*u.s,"redshift":_dl,"ngas":_dl,"nstar":_dl,"ndark":_dl,
         "Mgas":10**6*u.Msun,"Mstar":10**6*u.Msun,"Mdark":10**6*u.Msun, "MHI":10**6*u.Msun,
         "r_e_L":u.kpc,"r_e_M":u.kpc,"r_e_M_DM":u.kpc,
         "velDisp_star":u.km/u.s, "velDisp_DM":u.km/u.s,
         "kinE_star":u.km/u.s,"kinE_DM":10**6*u.Msun*u.km/u.s,
         "SFR":u.Msun/u.yr,
         "M_U":u.mag,"M_B":u.mag,"M_V":u.mag,"M_R":u.mag,"M_I":u.mag,"M_J":u.mag,"M_H":u.mag,"M_K":u.mag,
         "metall_star":_dl,"FeH_star":_dl,"MgFe_star":_dl,"metall_lum":_dl,"FeH_lum":_dl,"age_lum":_dl,
         "Mgas_re":10**6*u.Msun,"Mstar_re":10**6*u.Msun,"Mdark_re":10**6*u.Msun,
         "rcom_x":u.km,"rcom_y":u.km,"rcom_z":u.km,"vcom_x":u.km/u.s,"vcom_y":u.km/u.s,"vcom_z":u.km/u.s,
         "L_x":10**6*u.Msun*u.kpc*u.km/u.s,"L_y":10**6*u.Msun*u.kpc*u.km/u.s,"L_z":10**6*u.Msun*u.kpc*u.km/u.s,
         "r0":u.kpc,"sig(r0)":u.kpc,"n":_dl,"sig(n)":_dl,"mu0":u.mag*u.arcsec**2,"sig(mu0)":u.mag*u.arcsec**2}

meta = {"metall_star":"Metallicity weighted with star mass",
        "metall_lum":"Metallicity weighted with star luminosity (typically higher than metall_star)"}

exclude_names = ["flatStar","flatDM","bound_star","bound_DM", "I0", 'sig(I0)']

def get_sumfile(sumfilepath, **kwargs):
    table = ascii.read(sumfilepath, format='commented_header', exclude_names=exclude_names, **kwargs)
    for column, unit in units.items():
        table[column].unit = unit
    table.meta = meta
    return table

if __name__ == '__main__':
    datpath = "/home/michele/sim/MoRIA/results/sumfiles/69002.dat"
    table = get_sumfile(datpath)