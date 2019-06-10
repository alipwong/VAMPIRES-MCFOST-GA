import gaga as ga
import VAMPIRES_data
import star
import numpy as np

gene_definition = {"dust_mass": (1e-8, 3e-5),
                   "inner_radius": (1, 10),
                   "radius": (1, 800)}

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

    free_params = {"dust_mass": dm, "Rin": Rin, "radius":r}
    s = star.Star(free_params, default_params, WAVELENGTH, MCFOST_path, VAMP_data, load=False, verbose=False)
    s.calculate_reduced_chi2error()

    individual.data["data_Q"] = s.data_Q
    individual.data["data_U"] = s.data_U

    # We want to enforce the amplitude

    Q_obs_amp = max(VAMP_data.vhvv) - min(VAMP_data.vhvv) # around 0.1
    Q_sim_amp = max(s.data_Q) - min(s.data_Q) # previously around 0
    # general fitness is around 21 without it
    # really want to enforce amplitude

    llambda = 1000

    individual.fitness_score = s.reduced_chi2err + llambda * abs(Q_obs_amp - Q_sim_amp)

sim = ga.ga(gene_definition,
            evaluate,
            mutate = 0.6,
            sigma = 0.3,
            mutate_crossover = True,
            population_size = 25,
            verbose = True,
            results_folder = "amp_1000")

sim.run_simulation()

