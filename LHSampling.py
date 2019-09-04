import numpy as np
from smt.sampling_methods import LHS
import matplotlib.pyplot as plt
import time

from multistrain_model import simulate

# parameter range
CONTACTS_PER_HOST = [4, 12]
MU = [1/10, 1/3]  # recovery probability
SIGMA = [1/30, 1/10]  # immunity lost probability
BETA = [0.2, 0.8]  # infection probability
R = [0.01, 0.1]  # recombination probability per allele
TAO = [0.001, 0.005]  # mutation probability per allele
GAMMA = [0.02, 4]  # cross-immunity
N_LOCI = [2, 4]
N_NODES = [100, 500]

# fixed parameters
SEEDS_PER_STRAIN = 2
RANDOM = 1  # edge reconnecting probability of 1
REGULAR = 0     # edge reconnecting probability of 0
N_STEPS = 3000

# create multi-dimension parameter space
parameter_space = np.array([
    CONTACTS_PER_HOST,
    MU,
    SIGMA,
    BETA,
    R,
    TAO,
    GAMMA,
    N_LOCI,
    N_NODES
])

# LHSampling
sampling = LHS(xlimits=parameter_space)
n_sample_points = 1000
parameter_samples = sampling(n_sample_points)

# change numpy.array to python list
parameter_sample_list = parameter_samples.tolist()
# round int parameters to int
for sample in parameter_sample_list:
    sample[0] = int(round(sample[0]))
    sample[7] = int(round(sample[7]))
    sample[8] = int(round(sample[8]))

# diversity and discordance records
regular_diversity = []
regular_discordance = []
random_diversity = []
random_discordance = []

# simulate
i = 0
start_time = time.process_time()
for sample in parameter_sample_list:
    ran_div, ran_disc = simulate(*sample, SEEDS_PER_STRAIN, RANDOM, N_STEPS)
    reg_div, reg_disc = simulate(*sample, SEEDS_PER_STRAIN, REGULAR, N_STEPS)
    if ran_div and ran_disc and reg_div and reg_disc:
        print('---------------------', i, 'of', n_sample_points, '--------------------')
        i += 1
        random_diversity.append(ran_div)
        random_discordance.append(ran_disc)
        regular_diversity.append(reg_div)
        regular_discordance.append(reg_disc)
print('finished in', time.process_time() - start_time, 'seconds')

plt.subplot(121)
plt.scatter(random_diversity, regular_diversity)
plt.plot([0, 1], [0, 1])
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.xlabel('Diversity in random model')
plt.ylabel('Diversity in regular model')
plt.gca().set_aspect('equal', adjustable='box')

plt.subplot(122)
plt.scatter(random_discordance, regular_discordance)
plt.plot([0.4, 1], [0.4, 1])
plt.xlim(0.4, 1)
plt.ylim(0.4, 1)
plt.xlabel('Discordance in random model')
plt.ylabel('Discordance in regular model')
plt.gca().set_aspect('equal', adjustable='box')

plt.show()
