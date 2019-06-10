import numpy as np

def evaluate(individual):
	'''
		The Easom function has a minimim at (pi, pi), with a value of -1
	'''
	# constants


	#  unpack chromosome
	x = individual.CHROMOSOME['x']
	y = individual.CHROMOSOME['y']

	individual.FITNESS_SCORE = -np.cos(x) * np.cos(y) * np.exp(-((pow(x - np.pi, 2) + pow(y - np.pi, 2))))
