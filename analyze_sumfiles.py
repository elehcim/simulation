import os
import pandas as pd
import shutil
import matplotlib.pylab as plt
import numpy as np


header = ["# nr","time","redshift","ngas","nstar","ndark","Mgas","Mstar","Mdark","MHI","r_e_L","r_e_M","r_e_M_DM",
          "flatStar","flatDM","velDisp_star","velDisp_DM","kinE_star","kinE_DM","bound_star","bound_DM","SFR",
          "M_U","M_B","M_V","M_R","M_I","M_J","M_H","M_K","metall_star","FeH_star","MgFe_star","metall_lum","FeH_lum",
          "age_lum","Mgas_re","Mstar_re","Mdark_re","rcom_x","rcom_y","rcom_z","vcom_x","vcom_y","vcom_z",
          "L_x","L_y","L_z","I0","sig(I0)","r0","sig(r0)","n","sig(n)","mu0","sig(mu0)",
          "gas<1Re","gas<2Re","gas<3Re","gas<4Re","gas<5Re","gas<10Re","gas<30Re","gas all",
          "DM<1Re","DM<2Re","DM<3Re","DM<4Re","DM<5Re","DM<10Re","DM<30Re","DM all"]

def change_header(filename):
    backup = filename + ".orig"
    if os.path.isfile(backup):
        raise IOError("{} already exists".format(backup))
    shutil.copyfile(filename, backup)
    with open(backup, "r") as b:
        content = b.readlines()
        with open(filename, "w") as f:
            content[0] = content[0].replace("# nr", "   #").replace("gas all", "gas_all").replace("DM all", "DM_all")
            # print(content[0])
            for line in content:
                f.write(line)


def get_sumfile(sumfilepath):
    try:
        change_header(sumfilepath)
    except IOError:
        print("Header already modified")
    return pd.read_csv(sumfilepath, delim_whitespace=True, index_col=0)


# def ()
if __name__ == '__main__':
    datpath = "/home/michele/sim/MySimulations/Moria8Gyr_tidal/sim60003/60003.dat"
    df = get_sumfile(datpath)