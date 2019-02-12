# This script is used to experiment with feed a density file to MCFOST
# cmd: mcfost <parameter_file> ­density_file <your_density_file.fits.gz> ­3D (+ any other option)
# density dimensions: (1:n_rad, 1:nz, 1:n_az, 1:n_grains)

import sys
import numpy as np
import subprocess
import helper
import density_files
from astropy.io import fits
import time

mm_path = "../make-measurements-code/"   # path to where all the make-measurements code is
sys.path.insert(0, mm_path)
import star
import VAMPIRES_data

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

n_rad = 100
nz = 70
n_az = 1
n_grains = 1

default_params = "default_parameters_star"
free_params = {"image_symmetry":'F',
               "central_symmetry":'F',
               "axial_symmetry":'F',
               "Rout": 10,
               "n_incl": 2,
               "imax":90,
               "n_az":2,
               "n_az_angles":2,
               "az_max":90
               }

# -----
# CLEAR WORKSPACE
# -----

helper.clear_directory(mcfost_path)



# CREATE DENSITY FILE
# -----

# generate a random density file
# density_values = np.random.rand(n_rad, nz, n_az, n_grains)
# density_values = np.random.rand(1, 1, 70, 100) * 1000

# density file top and bottom
density_values = density_files.cross_section_bottom(nz, n_rad)
print(density_values)
print(density_values.shape)

hdu = fits.PrimaryHDU(density_values)
hdu.writeto(mcfost_path + density_file, overwrite=True)

# -----
# GENERATE STAR
# -----

s = star.Star(free_params, default_params, WAVELENGTH, mcfost_path, VAMP_data, load = False, verbose = False, density_file=density_file)
s.dm_data = s.load_data("data_disk/", "dust_mass_density.fits.gz")

