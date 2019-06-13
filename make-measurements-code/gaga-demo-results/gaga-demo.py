import gaga as ga
import VAMPIRES_data
import star
import numpy as np

gene_definition = {"dust_mass": (1e-7, 6e-7),
                   "inner_radius": (1, 20),
                   "radius": (1, 800),
                   "outer_radius": (1, 30)}

def evaluate(individual):

    WAVELENGTH = 0.75 # microns
    MCFOST_path = "../mcfost-workspace/"

    default_params = "default_parameters_star"

    VAMP_path = ""
    VAMP_data = "diffdata_RLeo_03_20170313_750-50_18holeNudged_0_0.idlvar"
    VAMP_data_info = "cubeinfoMar2017.idlvar"

    VAMP_data = VAMPIRES_data.VAMPIRES_data(VAMP_path + VAMP_data, VAMP_path + VAMP_data_info)

    # remove out bias from the observed data i.e. move the average to 1

    VAMP_data.vhvv -= np.mean(VAMP_data.vhvv) - 1
    VAMP_data.vhvvu -= np.mean(VAMP_data.vhvvu) - 1

    dm = individual.genes["dust_mass"]
    Rin = individual.genes["inner_radius"]
    r = individual.genes["radius"]
    Rout = individual.genes["outer_radius"]

    free_params = {"dust_mass": dm, "Rin": Rin, "radius":r, "Rout":Rout}
    s = star.Star(free_params, default_params, WAVELENGTH, MCFOST_path, VAMP_data, load=False, verbose=False)

    s.calculate_reduced_chi2error()

    individual.data["data_Q"] = s.data_Q
    individual.data["data_U"] = s.data_U
    individual.data["reduced_chi2err_Q"] = s.reduced_chi2err_Q
    individual.data["reduced_chi2err_U"] = s.reduced_chi2err_U
    individual.data["reduced_chi2err"] = s.reduced_chi2err

    # We want to enforce the standard deviation

    llambda = 10000 # -> if Q has a standard deviation of ~0, this should map to a fitness of ~150

    individual.fitness_score = s.reduced_chi2err

sim = ga.ga(gene_definition,
            evaluate,
            mutate = 0.6,
            sigma = 0.3,
            mutate_crossover = True,
            population_size = 25,
            verbose = True,
            results_folder = "init-vamp-fit")

sim.run_simulation()

