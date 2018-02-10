#! /usr/bin/env python3

import pandas as pd

header = ["time", "totIntEnergy", "totPotEnergy", "totKinEnergy", 
"IntEnergy_gas", "PotEnergy_gas", "KinEnergy_gas", 
"IntEnergy_dm", "PotEnergy_dm", "KinEnergy_dm", 
"IntEnergy_disk", "PotEnergy_disk", "KinEnergy_disk", 
"IntEnergy_bulge", "PotEnergy_bulge", "KinEnergy_bulge", 
"IntEnergy_s", "PotEnergy_s", "KinEnergy_s", 
"IntEnergy_bndry", "PotEnergy_bndry", "KinEnergy_bndry", 
"mass_gas", "mass_dm", "mass_disk", "mass_bulge", "mass_s", "mass_bndry"]

en = pd.read_csv("energy.txt", names=header, sep=' ')

# TODO add a function to parse the energy 