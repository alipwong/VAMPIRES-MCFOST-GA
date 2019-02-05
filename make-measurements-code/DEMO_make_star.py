import numpy as np
import star
import VAMPIRES_data
import matplotlib.pyplot as plt

free_params = {}

default_params = "default_parameters_star"

WAVELENGTH = 0.75	# microns
# location where the parameter file will be built and where the data folders will be created
MCFOST_path = "../mcfost-workspace/"

VAMP_path = ""
VAMP_data = "diffdata_RLeo_03_20170313_750-50_18holeNudged_0_0.idlvar"
VAMP_data_info = "cubeinfoMar2017.idlvar"

# load in external data
VAMP_data = VAMPIRES_data.VAMPIRES_data(VAMP_path + VAMP_data, VAMP_path + VAMP_data_info)
# remove out bias from the observed data i.e. move the average to 1

VAMP_data.vhvv -= np.mean(VAMP_data.vhvv) - 1
VAMP_data.vhvvu -= np.mean(VAMP_data.vhvvu) - 1

s = star.Star(free_params, default_params, WAVELENGTH, MCFOST_path, VAMP_data, load = False, verbose = True)

# s.display_data()

# plt.show()

# Additional code:
# s.sanity_check_star(shell = True)
# s.sanity_check_shell()
# s.display_PVRs()
# s.display_IQUV(interpolation = "nearest")
# print("Reduced chi^2 error:\t", s.reduced_chi2err)
# print(s.parameters.__dict__)
