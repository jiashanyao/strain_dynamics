import numpy as np
from smt.sampling_methods import LHS
import matplotlib.pyplot as plt
import time, multiprocessing, json

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

# result lists
ran_result_list = []
re_result_list = []


def log_ran_result(result):
    if result:
        ran_result_list.append(result)


def log_re_result(result):
    if result:
        re_result_list.append(result)


if __name__ == '__main__':
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
    n_sample_points = 5
    parameter_samples = sampling(n_sample_points)

    # change numpy.array to python list
    parameter_sample_list = parameter_samples.tolist()
    # round int parameters to int
    for sample in parameter_sample_list:
        sample[0] = int(round(sample[0]))
        sample[7] = int(round(sample[7]))
        sample[8] = int(round(sample[8]))

    # simulate
    i = 0
    start_time = time.process_time()
    pool = multiprocessing.Pool()
    for sample in parameter_sample_list:
        pool.apply_async(simulate, args=(*sample, SEEDS_PER_STRAIN, RANDOM, N_STEPS), callback=log_ran_result)
        pool.apply_async(simulate, args=(*sample, SEEDS_PER_STRAIN, REGULAR, N_STEPS), callback=log_re_result)
    pool.close()
    pool.join()
    print('finished in', time.process_time() - start_time, 'seconds')
    print(ran_result_list)
    print(re_result_list)
    with open('ran_result.json', 'w') as outfile:
        json.dump(ran_result_list, outfile, indent=2)
    with open('re_result.json', 'w') as outfile:
        json.dump(re_result_list, outfile, indent=2)

    # plt.subplot(121)
    # plt.scatter(random_diversity, regular_diversity)
    # plt.plot([0, 1], [0, 1])
    # plt.xlim(0, 1)
    # plt.ylim(0, 1)
    # plt.xlabel('Diversity in random model')
    # plt.ylabel('Diversity in regular model')
    # plt.gca().set_aspect('equal', adjustable='box')
    #
    # plt.subplot(122)
    # plt.scatter(random_discordance, regular_discordance)
    # plt.plot([0.4, 1], [0.4, 1])
    # plt.xlim(0.4, 1)
    # plt.ylim(0.4, 1)
    # plt.xlabel('Discordance in random model')
    # plt.ylabel('Discordance in regular model')
    # plt.gca().set_aspect('equal', adjustable='box')
    #
    # plt.show()
