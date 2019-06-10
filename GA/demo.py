#This is for any tests I want to make'

import Settings
import GA
import numpy as np

np.random.seed(1)
settings = Settings.Settings(n = 1)
sim = GA.GA(settings)

sim.run_simulation()

