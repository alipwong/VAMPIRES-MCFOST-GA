# This script is used to investigate the outputs (mostly the dimensions) from mcfost when running with the -output_density_grid flag.
# In particular grid.fits and dust_mass_density.fits

mm_path = "../make-measurements-code/"   # path to where all the make-measurements code is

import sys
import numpy as np
sys.path.insert(0, mm_path)
import star
import VAMPIRES_data

default_params = "default_parameters_star"

WAVELENGTH = 0.75	# microns
# location where the parameter file will be built and where the data folders will be created
MCFOST_path = "../mcfost-workspace/"

VAMP_path = mm_path
VAMP_data = "diffdata_RLeo_03_20170313_750-50_18holeNudged_0_0.idlvar"
VAMP_data_info = "cubeinfoMar2017.idlvar"

# load in external data
VAMP_data = VAMPIRES_data.VAMPIRES_data(VAMP_path + VAMP_data, VAMP_path + VAMP_data_info)
# remove out bias from the observed data i.e. move the average to 1

VAMP_data.vhvv -= np.mean(VAMP_data.vhvv) - 1
VAMP_data.vhvvu -= np.mean(VAMP_data.vhvvu) - 1

# STAR 1
free_params = {'n_az': 3, 'n_rad':20, 'nz': 70, 'n_rad_in':10}
s = star.Star(free_params, default_params, WAVELENGTH, MCFOST_path, VAMP_data, load = False)
s.grid_data = s.load_data("data_disk/", "grid.fits.gz")
s.dm_data = s.load_data("data_disk/", "dust_mass_density.fits.gz")
print(s.grid_data.shape)
print(s.dm_data.shape)

