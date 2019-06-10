# import matplotlib
# # matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import numpy as np
# import seaborn as sns
# import colours

plt.style.use('ggplot')

class Results:

    def __init__(self):

        self.trial_N = None
        self.results_folder = None
        self.gene_names = None

        self.epochs = 0

        # saves additional data that is defined in the evaluation function
        self.data = {"history": [],
                     "fitness": [],
                     "genes": [],
                     "diversity": [],
                     "average_diversity": []}
        # All are lists of lists. The inner list will contain the values, and the outer list will represent each epoch
        # history: the individuals at each epoch
        # genes: a dictionary with each gene as a key. The corresponding values are that of the genes at each epoch
        # fitness: the fitness scores of each individual at each epoch

    def gene_mins(self, gene):
        return [min(g[gene]) for g in self.data["genes"]]

    def gene_maxs(self, gene):
        return [max(g[gene]) for g in self.data["genes"]]

    def gene_median(self, gene):
        return [np.median(g[gene]) for g in self.data["genes"]]

    def gene_global_min(self, gene):
        return min(self.gene_mins(gene))

    def gene_global_max(self, gene):
        return max(self.gene_maxs(gene))

    def fitness_maxs(self):
        '''returns the maximum fitness in each generation'''
        return [max(f) for f in self.data["fitness"]]

    def fitness_mins(self):
        '''returns the minimum fitness in each generation'''
        return [min(f) for f in self.data["fitness"]]

    def fitness_global_min(self):
        return min(self.fitness_mins())

    def fitness_global_max(self):
        return max(self.fitness_maxs())

    def create_fitness_df(self):
        pass

    def __calculate_limits(self, gene, margin = 0, i = 0):
        center = self.gene_median(gene)[-1]

        gene_maxs = max(self.gene_maxs(gene)[i:])
        gene_mins = min(self.gene_mins(gene)[i:])

        diameter = max(abs(gene_maxs - center), abs(gene_mins - center))
        if not isinstance(margin, int):
            margin = margin[i]

        min_val = center - (1 + margin) * diameter
        max_val = center + (1 + margin) * diameter

        return min_val, max_val

    def plot_2d_init(self):
        pass

    def plot_2d_frame(self, i, start, ax1, ax2, x_gene, y_gene, colormap, log_scale, fmin, fmax):

        x = self.data["genes"][i + start][x_gene]
        y = self.data["genes"][i + start][y_gene]

        fitness_score = self.data["fitness"][i + start]
        if fmin is None:
            fitness_min = self.fitness_global_min()
        else:
            fitness_min = fmin

        if fmax is None:
            fitness_max = self.fitness_global_max()
        else:
            fitness_max = fmax

        if log_scale:
            fitness_score = np.log(fitness_score)
            fitness_min = np.log(fitness_min)
            fitness_max = np.log(fitness_max)

        sct_plt1 = ax1.scatter(x, y, c = fitness_score, cmap = colormap, vmin = fitness_min, vmax = fitness_max, s = 3, alpha = 0.5)


        sct_plt2 = ax2.scatter(x, y, c=fitness_score, cmap=colormap, vmin=fitness_min, vmax=fitness_max, s=3, alpha=0.5)

        # to remove later
        # ax1.scatter(1, 1, c='k', marker='*', s=6)
        # ax2.scatter(1, 1, c='k', marker='*', s=6)

        ax1.scatter(3, 2, c = 'k', marker = '*', s = 6)
        ax2.scatter(3, 2, c='k', marker='*', s=6)

        ax1.scatter(-2.805118, 3.131312, c='k', marker='*', s=6)
        ax2.scatter(-2.805118, 3.131312, c='k', marker='*', s=6)

        ax1.scatter(-3.779310, -3.283186, c='k', marker='*', s=6)
        ax2.scatter(-3.779310, -3.283186, c='k', marker='*', s=6)
        ax1.scatter(3.584428, -1.848126, c='k', marker='*', s=6)
        ax2.scatter(3.584428, -1.848126, c='k', marker='*', s=6)
        if i == 0 and start == 0:
            cb = plt.colorbar(sct_plt1)
            cb.solids.set_edgecolor("face")
            if log_scale:
                cb.set_label("Fitness score (log scale)")
            else:
                cb.set_label("Fitness score")

    def plot_2d(self, x_gene, y_gene,
                colormap = cm.rainbow,
                log_scale = True,
                scale = 0.1,
                fmin = False,
                fmax = False,
                fps = 30,
                inset = False):
        '''If reframe is True, the second plot will try to automatically zoom in during the animation.
        If reframe is false, the second image will be an automatic constant zoom showing 10% of the first image. You can change the zoom factor by changing scale'''

        fig = plt.figure(figsize = (6, 3))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
        fig.subplots_adjust(left=0.12, bottom=0.2, right=0.92, top=0.93, wspace=0.4, hspace=None)
        ax1.set_xlabel(x_gene)
        ax1.set_ylabel(y_gene)
        ax2.set_xlabel(x_gene)
        ax2.set_ylabel(y_gene)

        x_range = self.gene_global_max(x_gene) - self.gene_global_min(x_gene)
        y_range = self.gene_global_max(y_gene) - self.gene_global_min(y_gene)
        ax1.set_xlim(self.gene_global_min(x_gene) - 0.1 * x_range, self.gene_global_max(x_gene) + 0.1 * x_range)
        ax1.set_ylim(self.gene_global_min(y_gene) - 0.1 * y_range, self.gene_global_max(y_gene) + 0.1 * y_range)


        ax1.set_xlim(-5, 5)
        ax1.set_ylim(-5, 5)
        if inset:
            x_min = inset[0][0]
            x_max = inset[0][1]
            y_min = inset[1][0]
            y_max = inset[1][1]
        else:
            x_min, x_max = self.__calculate_limits(x_gene, margin = [scale - 1])
            y_min, y_max = self.__calculate_limits(y_gene, margin = [scale - 1])
        ax2.set_xlim(x_min, x_max)
        ax2.set_ylim(y_min, y_max)

        mark_inset(ax1, ax2, loc1 = 2, loc2 = 3, ec = 'k')

        max_frames = 400 # total number of frames per gif
        n_animations = int(np.ceil(self.epochs/max_frames))
        for i in range(n_animations):
            if i == n_animations - 1:
                n_frames = self.epochs - i * max_frames
            else:
                n_frames = max_frames
            print(n_frames)

            start = i * max_frames
            ani = animation.FuncAnimation(fig, self.plot_2d_frame, fargs = (start, ax1, ax2, x_gene, y_gene, colormap, log_scale, fmin, fmax), frames = n_frames, init_func = self.plot_2d_init, interval = 5, repeat = False)
            ani.save("{}{}_{}_{}_{}_{}_progression.gif".format(self.results_folder, self.trial_N, x_gene, y_gene, scale, start), writer='imagemagick', fps=fps)
        # plt.show()

    def plot_fitness(self):
        plt.figure()
        plt.title("Minimum fitness score")
        plt.xlabel("Epoch")
        plt.ylabel("fitness")
        plt.plot(self.fitness_mins())
        plt.savefig(fname = "{}{}_fitness".format(self.results_folder, self.trial_N))
        plt.show()

    def plot_diversity(self, gene = None):
        plt.figure()
        plt.title("Diversity")
        plt.plot(self.data["average_diversity"], label = "Average", linewidth = 0.5)
        for gene in self.gene_names:
            gene_div = [div[gene] for div in self.data["diversity"]]
            plt.plot(gene_div, label = gene, linewidth = 0.5)
        plt.legend()
        plt.savefig(fname="{}{}_diversity".format(self.results_folder, self.trial_N))
        plt.show()

    def plot_3d(self, x_gene, y_gene, z_gene):
        pass

    def print_best(self, minimise = True):
        best = False
        for epoch in self.data["history"]:
            for ind in epoch:

                if best == False:
                    best = ind
                else:
                    if minimise:
                        if ind.FITNESS_SCORE < best.FITNESS_SCORE:
                            best = ind
                    else:
                        if ind.FITNESS_SCORE > best.FITNESS_SCORE:
                            best = ind
        print(best.chromosome)
        print(best.FITNESS_SCORE)