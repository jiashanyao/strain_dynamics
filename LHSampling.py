import numpy as np
from smt.sampling_methods import LHS

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
RANDOMNESS = 1  # host contact network randomness, edge reconnecting probability
N_STEPS = 2000

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
n_sample_points = 50
parameter_samples = sampling(n_sample_points)

# change numpy.array to python list
parameter_sample_list = parameter_samples.tolist()
# round int parameters to int
for sample in parameter_sample_list:
    sample[0] = int(round(sample[0]))
    sample[7] = int(round(sample[7]))
    sample[8] = int(round(sample[8]))
# simulate
simulate(*parameter_sample_list[0], RANDOMNESS, N_STEPS)
