# This script is used to experiment with feed a density file to MCFOST
# cmd: mcfost <parameter_file> ­density_file <your_density_file.fits.gz> ­3D (+ any other option)
# density dimensions: (1:n_rad, 1:nz, 1:n_az, 1:n_grains)
import time
import sys
import numpy as np
import subprocess

import density_files
from astropy.io import fits
import time
import matplotlib.pyplot as plt

mm_path = "../make-measurements/"   # path to where all the make-measurements code is
sys.path.insert(0, mm_path)
import star
import helper
import VAMPIRES_data
start = time.time()
# -----
# DEFINE CONSTANTS AND LOAD EXTERNAL DATA
# -----
WAVELENGTH = 0.75   # in microns
mcfost_path = "../mcfost-workspace/"
density_file = "density.fits"

VAMP_path = mm_path
VAMP_data = "diffdata_RLeo_03_20170313_750-50_18holeNudged_0_0.idlvar"
VAMP_data_info = "cubeinfoMar2017.idlvar"

# load in external data
VAMP_data = VAMPIRES_data.VAMPIRES_data(VAMP_path + VAMP_data, VAMP_path + VAMP_data_info)
# remove out bias from the observed data i.e. move the average to 1

VAMP_data.vhvv -= np.mean(VAMP_data.vhvv) - 1
VAMP_data.vhvvu -= np.mean(VAMP_data.vhvvu) - 1

# -----
# DEFINE STAR
# -----
spot1 = (-1, 0.3, 45, 45) #amp, sigma, theta, phi (degrees)
spot2 = (30, 0.05, 180, 90)
# spot2 = (-1, 1, 30, 30)
default_params = "default_parameters_star"
free_params = {"nbr_photons_image":1.28e6,
               "dust_mass":1e-6,
               "Rin": 3,
               "Rout" : 30,
               "n_grains": 6,
               "amin":50,
               "amax":250,
               "n_az":200,
               "asp":1.5,
               "az_min":0,
               "az_max":0,
               "n_az_angles":1,
               "imin":90,
               "imax":90,
               "n_incl":1}#,
               #"spots":[spot1]}#, spot2]}

# -----
# CLEAR WORKSPACE
# -----

helper.clear_directory(mcfost_path)

# -----
# CREATE DENSITY FILE2
# -----


s = star.Star(free_params, default_params, WAVELENGTH, mcfost_path, VAMP_data, load = False, verbose = True)

end = time.time()
print("time")

print(end - start)

s.display_data(option = 1)
s.display_IQUV(scale = "off")
s.display_P(scale = "off")


plt.show()
