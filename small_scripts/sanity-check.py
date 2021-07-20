# This is a minimum viable product to demonstrate how to run MCFOST and obtain visabilities
import sys
import numpy as np
import matplotlib
matplotlib.use('TkAgg')

mm_path = "../make-measurements/"   # path to where all the make-measurements code is
sys.path.insert(0, mm_path)
import star
import star
import VAMPIRES_data
import helper
# -----
# DEFINE CONSTANTS AND LOAD EXTERNAL DATA
# -----
WAVELENGTH = 0.75   # in microns
mcfost_path = "../mcfost-workspace/"
density_file = "density.fits"
default_params = "sanity_check_parameters"
data_folder = "../temp/"

VAMP_path = mm_path
VAMP_data = "diffdata_RLeo_03_20170313_750-50_18holeNudged_0_0.idlvar"
VAMP_data_info = "cubeinfoMar2017.idlvar"

# load in external data
VAMP_data = VAMPIRES_data.VAMPIRES_data(VAMP_path + VAMP_data, VAMP_path + VAMP_data_info)
# remove out bias from the observed data i.e. move the average to 1

VAMP_data.vhvv -= np.mean(VAMP_data.vhvv) - 1
VAMP_data.vhvvu -= np.mean(VAMP_data.vhvvu) - 1

# -----
# DEFINE PARAMETERS
# -----

free_params = {"Rin": 4, "Rout": 4.1, "dust_mass": 1e-10, "radius": 800, "size": 501}

# helper.clear_directory(mcfost_path)
s = star.Star(free_params, default_params, WAVELENGTH, mcfost_path, VAMP_data, load = False, verbose = True)
s.sanity_check_star()
