import sys
import numpy as np
# from tabulate import tabulate
import Individual
import Results
import pickle
import bz2

class GA:

    def __init__(self, settings, load = False):

        # Import the evaluation function used to calculate the fitness score
        self.EVAL_FILE = __import__(settings.EVAL_FILENAME)
        self.EVALUATE = self.EVAL_FILE.evaluate

        # Gene constants

        self.NUMBER_GENES = len(settings.GENE_NAMES)
        self.GENE_NAMES = settings.GENE_NAMES
        self.LOWER_BOUNDS = {settings.GENE_NAMES[i]: settings.LOWER_BOUNDS[i] for i in range(self.NUMBER_GENES)}
        self.UPPER_BOUNDS = {settings.GENE_NAMES[i]: settings.UPPER_BOUNDS[i] for i in range(self.NUMBER_GENES)}
        self.LOGSPACED = settings.LOGSPACED

        # Simulation constants
        self.PROBLEM_TYPE = settings.PROBLEM_TYPE
        self.REVERSE = self.PROBLEM_TYPE == "maximise"
        self.POPULATION_SIZE = settings.POPULATION_SIZE

        self.cloneN = settings.clone     # percentage of the new generation produced by cloning
        self.mutateN = settings.mutant   # percentage of the new generation produced by mutation
        self.crossoverN = self.POPULATION_SIZE - self.cloneN - self.mutateN # number produced by crossover

        if self.crossoverN < 0:
            sys.exit("""...
ERROR in settings.
The clone and mutant rates are too high. No individuals produced by crossover.
Reduce either the clone or mutant rate (or both).
...""")

        self.genePoolN = 2 * self.crossoverN

        # Functions
        self.validate_chromosome = settings.validate_chromosome

        if settings.SELECTION == "tournament":
            self.SELECTION = self.tournament
        if settings.SELECTION == "rouletteWheel":
            self.SELECTION = self.rouletteWheel

        if settings.CROSSOVER == "multiUniformBrood":
            self.CROSSOVER = self.multiUniformBrood
        if settings.CROSSOVER == "multiUniform":
            self.CROSSOVER = self.multiUniform
        if settings.CROSSOVER == "arithmetic":
            self.CROSSOVER = self.arithmetic

        if settings.MUTATION == "random":
            self.MUTATION = self.randomMutation

        if settings.MUTATION == "gaussian":
            self.MUTATION = self.GaussianMutation
            self.SIGMA = settings.SIGMA

        # Assign constants to the Individual class
        Individual.Individual.EVALUATE = self.EVALUATE
        Individual.Individual.VALIDATE = self.validate_chromosome

        # Hyper parameters
        self.winRate = settings.winRate

        # termination condition
        self.EPOCH = settings.EPOCH

        # prints out population details at each epoch
        self.EXPLICIT = settings.EXPLICIT

        self.mutate_crossover = settings.mutate_crossover
        self.cross_scale = settings.cross_scale
        self.speciate = settings.speciate

        self.RESULTS_FOLDER = settings.RESULTS_FOLDER
        self.RESULTS = Results.Results()
        self.RESULTS.gene_names = settings.GENE_NAMES
        self.RESULTS.results_folder = settings.RESULTS_FOLDER
        self.RESULTS.trial_N = settings.TRIAL_N

        # Name of the folder of the simulation we want to continue from
        # Otherwise False
        self.HISTORY_FOLDER = load


        # Keep a record of the previous simulations
        if load == False:
            with open(self.RESULTS_FOLDER + "history.txt", 'w') as f:
                f.write(self.RESULTS_FOLDER)
        else:
            with open(self.HISTORY_FOLDER + "history.txt", 'r') as f:
                hist = f.read()
            with open(self.RESULTS_FOLDER + "history.txt", 'w') as f:
                f.write(hist)
                f.write(self.RESULTS_FOLDER)

        with bz2.BZ2File(self.RESULTS_FOLDER + "settings_obj", 'wb') as f:
            pickle.dump(settings, f)

    def run_simulation(self):

        ''' Forced to actually run a simulation'''

        self.init_population()
        for i in range(self.EPOCH):

            if self.EXPLICIT:
                self.print_ind(self.newPopulation, title="EPOCH: " + str(i))
            self.save_results(i)

            self.REFRESH()
            self.SELECTION()
            if self.speciate:
                self.order_gene_pool(i)
            # Create a new population and populate with crossover, mutation and cloning
            self.newPopulation = []
            self.CROSSOVER()
            # self.print_ind(self.newPopulation, title="New population (after crossover)")
            self.CLONE()
            # self.print_ind(self.newPopulation, title="New population (after cloning)")
            self.MUTATION()
        if self.EXPLICIT:
            self.print_ind(self.newPopulation, title="final population")
        self.save_results(self.EPOCH)


    def measure_diversity(self):
         # coefficient of variation (relative standard deviation)
        div = {}
        for gene in self.GENE_NAMES:
            # all the values of a particular gene in the population
            vals = [ind.CHROMOSOME[gene] for ind in self.newPopulation]
            div[gene] = np.std(vals)/abs(np.mean(vals))

        self.diversity = div
        self.average_diversity = np.mean(list(div.values()))

    def tune_hyper_parameters(self):
        ''' Tune the hyper-parameters based on the diversity of the population'''
        pass

    def random_sample_logspace(self, start, end):
        sample = np.random.uniform(np.log(start), np.log(end))
        return np.exp(sample)

    def init_population(self):

        ''' Create an initial population within the specified bounds '''

        # Start with a random population
        if self.HISTORY_FOLDER == False:
            # Create a random set of chromosomes
            init_chromosomes = []
            for individual in range(self.POPULATION_SIZE):
                # create a random chromosome
                chromosome = {}
                for gene in self.GENE_NAMES:
                    if gene in self.LOGSPACED:
                        chromosome[gene] = self.random_sample_logspace(self.LOWER_BOUNDS[gene], self.UPPER_BOUNDS[gene])
                    else:
                        chromosome[gene] = np.random.uniform(self.LOWER_BOUNDS[gene], self.UPPER_BOUNDS[gene])
                init_chromosomes.append(chromosome)

            # Create individuals
            self.newPopulation = [Individual.Individual(chromosome) for chromosome in init_chromosomes]
        else:
            with bz2.BZ2File(self.HISTORY_FOLDER + "results_obj", 'rb') as f:
                results = pickle.load(f)
                self.newPopulation = results.data["history"][-1]

    def REFRESH(self):
        ''' Refreshes the world after each epoch'''
        self.population = self.newPopulation
        self.newPopulation = []

#   SELECTION
#   Each selection method creates self.genePool through selection from self.population which
#   is a list of individuals that will go on to reproduce. The number of individuals is
#   specified by self.genePoolN

    def tournament(self):
        ''' Tournament selection. Two individuals are selected to enter into the tournament.
		The individual with the higher fitness score has a winRate chance of winning the tournament, while the individual with the lower fitness score has a (1 - winRate) chance of winning tournament. The tournament winner goes into the gene pool
		'''

        self.genePool = []  # create an empty gene pool

        # 	each tournament adds one individual to the gene pool
        for tournament in range(self.genePoolN):
            participants = np.random.choice(self.population, size = 2, replace = False)

            # sort participants by rank
            participants = sorted(participants, key = lambda individual: individual.FITNESS_SCORE, reverse = self.REVERSE)


            if np.random.uniform() < self.winRate:
                winner = participants[0]
            else:
                winner = participants[1]

            self.genePool.append(winner)

    def rouletteWheel(self):

        max_fitness = max([ind.FITNESS_SCORE for ind in self.population])

        inv_fitness = [max_fitness -  ind.FITNESS_SCORE + 1 for ind in self.population]

        cd = [sum(inv_fitness[:(i+1)]) for i in range(len(inv_fitness))]
        cd = np.asarray(cd)

        ind_selection = [sum(cd < cd[-1]*np.random.random()) for i in range(self.genePoolN)]

        self.genePool = [self.population[i] for i in ind_selection]

#   CROSSOVER
#   Each crossover method uses individuals from self.genePool and adds self.crossoverN individuals to self.newPopulation.
#   Note: self.genePoolN is always even

    def order_gene_pool(self, epoch):

        # select gene:
        gene_order = self.GENE_NAMES[epoch % len(self.GENE_NAMES)]
        self.genePool.sort(key = lambda ind: ind.CHROMOSOME[gene_order])

    def arithmetic(self):
        r = np.random.uniform()
        for pair in range(0, 2*self.crossoverN, 2):
            mother = self.genePool[pair]
            father = self.genePool[pair + 1]

            chromosome = {}
            for gene in self.GENE_NAMES:
                chromosome[gene] = r * mother.CHROMOSOME[gene] + (1-r) * father.CHROMOSOME[gene]

            # create individuals
            child = Individual.Individual(chromosome)

            # take the best child and add it to the next generation
            self.newPopulation.append(child)

    def multiUniform(self):
        '''This creates 1 child using the mask. The mask created is uniformly and randomly chosen.'''

        offspring = []

        for pair in range(0, 2*self.crossoverN, 2):
            mother = self.genePool[pair]
            father = self.genePool[pair + 1]

            # size: no. of parameters it picks. At least 1 gene, but less than self.NUMBER_GENES parameters
            # chooses randomly without replacement which parameters to put in the mask
            mask = np.random.choice(self.GENE_NAMES, size = np.random.randint(1, self.NUMBER_GENES), replace = False)

            chromosome = {}
            for gene in self.GENE_NAMES:
                if gene in mask:
                    chromosome[gene] = mother.CHROMOSOME[gene]
                else:
                    chromosome[gene] = father.CHROMOSOME[gene]

            if self.mutate_crossover:
                g = np.random.choice(self.GENE_NAMES)
                chromosome[g] = np.random.normal(loc = chromosome[g], scale = abs(self.cross_scale * chromosome[g]))

            # create individuals
            child = Individual.Individual(chromosome)

            # take the best child and add it to the next generation
            self.newPopulation.append(child)

    def multiUniformBrood(self):
        '''This creates 2 children using the mask and the inverse mask. We then take the best child. The mask created is uniformly and randomly chosen.'''

        offspring = []

        for pair in range(0, 2*self.crossoverN, 2):
            mother = self.genePool[pair]
            father = self.genePool[pair + 1]

            # size: no. of parameters it picks. At least 1 gene, but less than self.NUMBER_GENES parameters
            # chooses randomly without replacement which parameters to put in the mask
            mask = np.random.choice(self.GENE_NAMES, size = np.random.randint(1, self.NUMBER_GENES), replace = False)

            chromosome_1 = {}
            chromosome_2 = {}
            for gene in self.GENE_NAMES:
                if gene in mask:
                    chromosome_1[gene] = mother.CHROMOSOME[gene]
                    chromosome_2[gene] = father.CHROMOSOME[gene]
                else:
                    chromosome_1[gene] = father.CHROMOSOME[gene]
                    chromosome_2[gene] = mother.CHROMOSOME[gene]

            if self.mutate_crossover:
                for gene in self.GENE_NAMES:
                    chromosome_1[gene] = np.random.normal(loc = chromosome_1[gene], scale = self.cross_scale * chromosome_1[gene])
                    chromosome_2[gene] = np.random.normal(loc=chromosome_2[gene], scale=self.cross_scale * chromosome_2[gene])

            if self.mutate_crossover:
                g1 = np.random.choice(self.GENE_NAMES)
                chromosome_1[g1] = np.random.normal(loc = chromosome_1[g1], scale = abs(self.cross_scale * chromosome_1[g1]))
                g2 = np.random.choice(self.GENE_NAMES)
                chromosome_2[g2] = np.random.normal(loc=chromosome_2[g2], scale=abs(self.cross_scale * chromosome_2[g2]))

            # create individuals
            child_1 = Individual.Individual(chromosome_1)
            child_2 = Individual.Individual(chromosome_2)

            children = sorted([child_1, child_2], key = lambda individual: individual.FITNESS_SCORE, reverse = self.REVERSE)

            # take the best child and add it to the next generation
            self.newPopulation.append(children[0])

#   MUTATE
    def randomMutation(self):
        ''' This is the random mutation function described in Tang, Tseng, 2012 '''

        # Select individuals from the previous population to mutate. an individual can only be selected once
        to_mutate = np.random.choice(self.population, size = self.mutateN, replace = False)

        # Work out the maximum and minimum value of each parameter
        # A matrix where each row are the genes of an individual
        gene_max = {gene: max([i.CHROMOSOME[gene] for i in self.population]) for gene in self.GENE_NAMES}
        gene_min = {gene: min([i.CHROMOSOME[gene] for i in self.population]) for gene in self.GENE_NAMES}

        for ind in to_mutate:
            chr = ind.CHROMOSOME # original chromosome
            delta = {gene: max(2 * (chr[gene] - gene_min[gene]), 2 * (gene_max[gene] - chr[gene])) for gene in self.GENE_NAMES}
            new_chr = {gene: chr[gene] + delta[gene] * (np.random.uniform(0,1) - 0.5) for gene in self.GENE_NAMES}
            mutant = Individual.Individual(new_chr)
            self.newPopulation.append(mutant)

    def GaussianMutation(self, replace = False):

        to_mutate = np.random.choice(self.population, size=self.mutateN, replace=False)

        for ind in to_mutate:
            chr = ind.CHROMOSOME
            chr_copy = {k:v for k, v in chr.items()}

            # pick a gene to mutate
            g = np.random.choice(self.GENE_NAMES)
            chr_copy[g] = min(max(np.random.normal(loc = chr_copy[g], scale = self.SIGMA), self.LOWER_BOUNDS[g]), self.UPPER_BOUNDS[g])
            mutant = Individual.Individual(chr_copy)

            self.newPopulation.append(mutant)


#   CLONE
    def CLONE(self):

        # rank the previous population
        ranked = sorted(self.population, key = lambda individual: individual.FITNESS_SCORE, reverse = self.REVERSE)

        for i in range(self.cloneN):
            self.newPopulation.append(ranked[i])

    # PRINTING

    def print_ind(self, list_ind, dp = ".4f", title = None):
        ''' Convenience function. Pass in a list of individuals and it will print out the chromosomes and the fitness score'''
        if title:
            print("=" * len(title))
            print(title.upper())

        data = [[i.FITNESS_SCORE] + [i.CHROMOSOME[gene] for gene in self.GENE_NAMES] for i in list_ind]

        # print(tabulate(data, headers = ["Fitness Score"] + self.GENE_NAMES, tablefmt = "rst", floatfmt = dp))
        for i in data:
            print(i)

    def save_results(self, epoch):
        '''save results at each epoch'''
        self.measure_diversity()

        self.RESULTS.data["history"].append(self.newPopulation)
        self.RESULTS.data["fitness"].append([i.FITNESS_SCORE for i in self.newPopulation])
        self.RESULTS.data["genes"].append({gene: [i.CHROMOSOME[gene] for i in self.newPopulation] for gene in self.GENE_NAMES})
        self.RESULTS.data["diversity"].append(self.diversity)
        self.RESULTS.data["average_diversity"].append(self.average_diversity)

        # data saved at the end of the evaluation function
        eval_data = self.newPopulation[0].data.keys()
        for i in eval_data:
            if i not in self.RESULTS.data:
                self.RESULTS.data[i] = []

        for j in eval_data:
            self.RESULTS.data[j].append([i.data[j] for i in self.newPopulation])

        self.RESULTS.epochs += 1

        with bz2.BZ2File(self.RESULTS_FOLDER + "results_obj", 'wb') as f:
            pickle.dump(self.RESULTS, f)