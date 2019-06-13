import bz2
import pickle
import numpy as np
import VAMPIRES_data
import matplotlib.pyplot as plt

def display_observed_data(VAMP_data, ax, max_baseline=8, Stokes='Q', ylims=False, marker='x'):
    ''' Displays the data from the observation. Specify whether you want the data from Stokes Q or Stokes U. The maximum baseline is the largest baseline that will be displayed.'''

    baselines = VAMP_data.blengths
    # only select the baselines that are smaller than the maximum baseline
    baseline_colours = baselines[baselines <= max_baseline]

    if Stokes == 'Q':
        data = VAMP_data.vhvv
        errors = VAMP_data.vhvverr
        title = "Observed Data (Stokes Q)"
    else:
        Stokes = 'U'
        data = VAMP_data.vhvvu
        errors = VAMP_data.vhvvuerr
        title = "Observed Data (Stokes U)"

    # Plots the markers
    scatter_plot = ax.scatter(VAMP_data.bazims, data, s=25, marker=marker, label="observed", c=baseline_colours)
    clb = plt.colorbar(scatter_plot, label="Baseline length (m)")
    bar_colour = clb.to_rgba(baseline_colours)

    if ylims:
        ax.set_ylim(ylims)

    # Plots the errorbars
    ax.errorbar(VAMP_data.bazims, data, yerr=errors, marker='', linestyle='', alpha=0.8, capsize=0, zorder=0,
                ecolor=bar_colour)

    ax.set_title(title)
    ax.set_xlabel("Baseline Azimuth Angle (rads)")
    ax.set_ylabel("Polarised Visibility Ratio")


def display_sim_data(VAMP_data, data, ax, max_baseline=8, Stokes='Q', ylims=False, colorbar=True, marker='x'):
    ''' Displays the simulated data produced from the model. Specify whether you want the data from Stokes Q or Stokes U. The maximum baseline is the largest baseline that will be displayed.'''

    baselines = VAMP_data.blengths
    # only select the baselines that are smaller than the maximum baseline
    baseline_colours = baselines[baselines <= max_baseline]

    if Stokes == 'Q':
        title = "Simulated Data (Stokes Q)"
    else:
        Stokes = 'U'
        title = "Simulated Data (Stokes U)"

    scatter_plot = ax.scatter(VAMP_data.bazims, data, s=25, marker=marker, label="simulated",
                              c=baseline_colours)

    if ylims:
        ax.set_ylim(ylims)


VAMP_path = ""
VAMP_data = "diffdata_RLeo_03_20170313_750-50_18holeNudged_0_0.idlvar"
VAMP_data_info = "cubeinfoMar2017.idlvar"

# load in external data
VAMP_data = VAMPIRES_data.VAMPIRES_data(VAMP_path + VAMP_data, VAMP_path + VAMP_data_info)
VAMP_data.vhvv -= np.mean(VAMP_data.vhvv) - 1
VAMP_data.vhvvu -= np.mean(VAMP_data.vhvvu) - 1

with bz2.open("./gaga-demos/fit_Q_new_bounds_100/results_obj", "rb") as f:
    results = pickle.load(f)

results.results_folder = "./gaga-demos/fit_Q_new_bounds_100"


# results.plot_fitness()
# results.animate("dust_mass", "inner_radius", fmax = 150)
# results.animate("dust_mass", "radius", fmax = 150)
# results.animate("radius", "inner_radius", fmax = 150)

print(len(results.data['data_Q'][0]))

ind = np.argmin(results.data['fitness'], axis = 1)

#
for i in [10]:
    Q = results.data['data_Q'][i][ind[i]]
    U = results.data['data_U'][i][ind[i]]
    fig = plt.figure(figsize=(12, 5))

    ax = fig.add_subplot(121)
    display_observed_data(VAMP_data, ax, Stokes = 'Q')
    display_sim_data(VAMP_data, Q, ax, Stokes="Q", marker = 'o')
    ax = fig.add_subplot(122)
    display_observed_data(VAMP_data, ax, Stokes = 'U')
    display_sim_data(VAMP_data, U, ax, Stokes="U", marker = 'o')
plt.savefig(results.results_folder + 'final_fit')
plt.show()
#
# print(np.std(Q))
# print(np.std(VAMP_data.vhvv))
# #
#
# print(VAMP_data.vhvv)
# print(max(VAMP_data.vhvv) - min(VAMP_data.vhvv))
# print(max(Q) - min(Q))