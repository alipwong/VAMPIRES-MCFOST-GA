import os

class Settings:

    def __init__(self, n = 0):

        # CREATE RESULTS FOLDER
        self.TRIAL_N = n
        self.PARENT_FOLDER = ""
        self.RESULTS_FOLDER = self.PARENT_FOLDER + "demo" + "/"
        self.TRIAL_DESC = "Modifying code to allow more than just the chromosome of each individual to be saved"

        if os.path.exists(self.RESULTS_FOLDER):
            input("""...
Are you sure you want to continue?
You will be overwriting data.
...""")
        else:
            os.makedirs(self.RESULTS_FOLDER)

        with open(self.RESULTS_FOLDER + "notes.txt", 'w') as f:
            f.write(self.TRIAL_DESC)

        self.EVAL_FILENAME = "Himmelblau_eval"  # name of the file containing the evaluation function

        self.GENE_NAMES = ["x", "y"]  # a list of all the parameters
        self.LOWER_BOUNDS = [-5, -5]  # must be in the same order as in gene_names
        self.UPPER_BOUNDS = [5, 5]
        self.LOGSPACED = []

        self.PROBLEM_TYPE_OPT = ["maximise", "minimise"]  # problem type
        self.PROBLEM_TYPE = self.PROBLEM_TYPE_OPT[1]

        # methods
        self.SELECTION_OPT = ["tournament", "rouletteWheel"]
        self.SELECTION = self.SELECTION_OPT[0]

        self.CROSSOVER_OPT = ["multiUniform", "multiUniformBrood", "arithmetic"]
        self.CROSSOVER = self.CROSSOVER_OPT[0]
        self.mutate_crossover = True   # if true, will mutate a single gene a small amount
        self.cross_scale = 0.01 # only used with mutate_crossover
        self.speciate = True

        self.MUTATION_OPT = ["random", "gaussian"]
        self.MUTATION = self.MUTATION_OPT[1]
        self.SIGMA = 0.5 # Only used for gaussian mutati

        # hyper parameters
        self.winRate = 0.7

        self.POPULATION_SIZE = 100

        # these values will be rounded up. To ensure the population size remains constant, the remaining individuals will be produced via crossover
        self.clone = 1      # Number cloned
        self.mutant = 40  # Number mutated

        # termination condition
        self.EPOCH = 50

        self.EXPLICIT = False    # prints out population details at each iteration



    def validate_chromosome(self, chromosome):
        ''' Ensure that a chromosome is valid. If invalid, returns a new valid chromosome. '''

        # keep chromosome within initial bounds
        for i in range(len(self.GENE_NAMES)):
            gene = self.GENE_NAMES[i]
            lb = self.LOWER_BOUNDS[i]
            ub = self.UPPER_BOUNDS[i]
            if chromosome[gene] < lb:
                chromosome[gene] = lb
            if chromosome[gene] > ub:
                chromosome[gene] = ub
        return chromosome
