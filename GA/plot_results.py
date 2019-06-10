import pickle
import bz2

RESULTS_FOLDER = "demo/"

with bz2.BZ2File(RESULTS_FOLDER + "results_obj", 'rb') as f:
    results = pickle.load(f)

results.plot_fitness()
# results.plot_diversity()
# results.plot_2d('x', 'y', scale = 0.1, log_scale = False, fmin = 0, fmax = 20, fps = 50, inset = [[0.9, 1.1], [0.9, 1.1]])
# results.plot_2d('x', 'y', scale = 0.1, log_scale = False, fmin = 0, fmax = 20, fps = 50, inset = [[-4, -3], [-4, -2.5]])
# results.plot_2d('x', 'y', scale = 0.1, log_scale = False, fmin = 0, fmax = 20, fps = 50, inset = [[-3, -2.5], [2.5, 4]])
results.plot_2d('x', 'y', scale = 0.1, log_scale = False, fmin = 0, fmax = 20, fps = 50)
print(results.data)

results.print_best()

