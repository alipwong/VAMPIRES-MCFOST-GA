import gaga as ga
import VAMPIRES_data
import star
import numpy as np

gene_definition = {"dust_mass": (1e-7, 6e-7),
                   "inner_radius": (1, 20),
                   "radius": (1, 800),
                   "outer_radius": (1, 30)}

def evaluate(individual):

    WAVELENGTH = 0.75 # wavelength of observation in microns

    # location where intermediate mcfost files are created
    MCFOST_path = "../mcfost-workspace/"

    # name of the python file containing the Default_parameters class
    default_params = "default_parameters_star"

    # path to the VAMPIRES data
    VAMP_path = ""
    VAMP_data = "diffdata_RLeo_03_20170313_750-50_18holeNudged_0_0.idlvar"
    VAMP_data_info = "cubeinfoMar2017.idlvar"

    VAMP_data = VAMPIRES_data.VAMPIRES_data(VAMP_path + VAMP_data, VAMP_path + VAMP_data_info)

    # remove bias from the observed data i.e. move the average to 1
    VAMP_data.vhvv -= np.mean(VAMP_data.vhvv) - 1
    VAMP_data.vhvvu -= np.mean(VAMP_data.vhvvu) - 1

    # unpack the genes and use them to define the parameters for the model
    dm = individual.genes["dust_mass"]
    Rin = individual.genes["inner_radius"]
    r = individual.genes["radius"]
    Rout = individual.genes["outer_radius"]

    free_params = {"dust_mass": dm, "Rin": Rin, "radius":r, "Rout":Rout}

    # run mcfost
    s = star.Star(free_params, default_params, WAVELENGTH, MCFOST_path, VAMP_data, load=False, verbose=False)

    # calculate the reduced chi^2 error and use it as the fitness score
    s.calculate_reduced_chi2error()
    individual.fitness_score = s.reduced_chi2err

    # Optional: storing some useful data
    individual.data["data_Q"] = s.data_Q
    individual.data["data_U"] = s.data_U
    individual.data["reduced_chi2err_Q"] = s.reduced_chi2err_Q
    individual.data["reduced_chi2err_U"] = s.reduced_chi2err_U
    individual.data["reduced_chi2err"] = s.reduced_chi2err

# NOTE: This GA has not been tuned
sim = ga.ga(gene_definition,
            evaluate,
            mutate = 0.6,
            sigma = 0.3,
            mutate_crossover = True,
            population_size = 25,
            verbose = True,
            results_folder = "gaga-demo-results")

sim.run_simulation()

