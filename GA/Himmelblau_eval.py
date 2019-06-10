import numpy as np
def evaluate(individual):

		#The Eggholder function has a minimum at (512, 404.2319)


	#  unpack chromosome
	x = individual.CHROMOSOME['x']
	y = individual.CHROMOSOME['y']

	individual.FITNESS_SCORE = pow(pow(x, 2) + y - 11, 2) + pow(x + pow(y, 2) - 7, 2)


