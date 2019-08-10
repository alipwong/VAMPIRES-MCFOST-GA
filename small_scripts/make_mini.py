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
import pickle
import itertools

mm_path = "../make-measurements/"   # path to where all the make-measurements code is
sys.path.insert(0, mm_path)
import star
import VAMPIRES_data

# -----
# DEFINE CONSTANTS AND LOAD EXTERNAL DATA
# -----
WAVELENGTH = 0.75   # in microns
mcfost_path = "../mcfost-workspace/"
density_file = "density.fits"
default_params = "default_parameters_star"
data_folder = "../new-sym/"

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

# start, end number of values
radius = np.geomspace(100, 1000, 10) # stellar radius
rin = np.geomspace(1, 10, 10)
rfrac = np.geomspace(1.5, 5, 10)
density = np.geomspace(0.05, 1, 10)


# -----
# PERFORM SWEEP
# -----

counter = 1
c = 0
params = list(itertools.product(radius, rin, rfrac, density))
n = len(params)

for radius, rin, rfrac, density in params:

    c += 1
    rout = rfrac * rin
    volume = (4/3 * np.pi * (pow(rout, 3) - pow(rin, 3)))
    dm = density * volume / 1e10

    print("{} of {}".format(c, n))

    helper.clear_directory(mcfost_path)

    print(dm, rin, rout, radius)
    free_params = {"dust_mass": dm, "Rin": rin, "radius": radius, "Rout":rout}

    s = star.Star(free_params, default_params, WAVELENGTH, mcfost_path, VAMP_data, load = False, verbose = False)
    if s.made:
        s.observe()
        vals = [free_params, s.data_Q, s.data_U]

        with open(data_folder + "{}.pkl".format(counter), 'wb') as f:
            pickle.dump(vals, f)

        counter += 1

# import time
#
# helper.clear_directory(mcfost_path)
# start = time.time()
# free_params = {"dust_mass": 1e-6}
#
# s = star.Star(free_params, default_params, WAVELENGTH, mcfost_path, VAMP_data, load = False, verbose = True)
# s.observe()
# vals = [free_params, s.data_Q, s.data_U]
#
# with open(data_folder + "{}.pkl".format(counter), 'wb') as f:
#     pickle.dump(vals, f)
# end = time.time()
#
# print((end - start))
