import multiprocessing, json, numpy
from smt.sampling_methods import LHS
from multistrain_model import simulate

# fixed parameters
SEEDS_PER_STRAIN = 1
RANDOM = 1  # edge reconnecting probability of 1
REGULAR = 0  # edge reconnecting probability of 0
N_STEPS = 3000
N_SAMPLE_POINTS = 1000


def log_result(result, result_list):
    if result:
        result_list.append(result)


def run_lhs(contacts_per_host, mu, sigma, beta, r, tau, gamma, n_loci, n_nodes, ran_filename, re_filename):
    # create multi-dimension parameter space
    parameter_space = numpy.array([
        contacts_per_host,
        mu,
        sigma,
        beta,
        r,
        tau,
        gamma,
        n_loci,
        n_nodes
    ])

    # LHSampling
    sampling = LHS(xlimits=parameter_space)
    parameter_samples = sampling(N_SAMPLE_POINTS)

    # change numpy.array to python list
    parameter_sample_list = parameter_samples.tolist()
    # round int parameters to int
    for sample in parameter_sample_list:
        sample[0] = int(round(sample[0]))
        sample[1] = 1/sample[1]
        sample[2] = 1/sample[2]
        # sample[4] = 0   # disable recombination
        sample[7] = int(round(sample[7]))
        sample[8] = int(round(sample[8]))

    # simulate and write results to file
    pool = multiprocessing.Pool()
    ran_result_list = []
    re_result_list = []
    for sample in parameter_sample_list:
        pool.apply_async(simulate, args=(*sample, SEEDS_PER_STRAIN, RANDOM, N_STEPS), callback=lambda result: log_result(result, ran_result_list))
        pool.apply_async(simulate, args=(*sample, SEEDS_PER_STRAIN, REGULAR, N_STEPS), callback=lambda result: log_result(result, re_result_list))
    pool.close()
    pool.join()
    with open(ran_filename, 'w') as outfile:
        json.dump(ran_result_list, outfile, indent=2)
    with open(re_filename, 'w') as outfile:
        json.dump(re_result_list, outfile, indent=2)


if __name__ == '__main__':
    # parameter range
    contacts_per_host = [4, 12]
    mu = [3, 10]  # recovery probability
    sigma = [10, 30]  # immunity lost probability
    beta = [0.2, 0.8]  # infection probability
    r = [0.01, 0.1]  # recombination probability per allele
    tau = [0.001, 0.005]  # mutation probability per allele
    gamma = [0.02, 4]  # cross-immunity
    n_loci = [2, 4]
    n_nodes = [100, 500]
    run_lhs(contacts_per_host, mu, sigma, beta, r, tau, gamma, n_loci, n_nodes,
            'ran_new_recombination.json', 're_new_recombination.json')
