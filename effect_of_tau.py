import multiprocessing
from matplotlib import pyplot as plt
from multistrain_model import simulate

min_tau = 0.001
max_tau = 0.005
n_samples = 10


def log_result(result, result_list):
    if result:
        result_list.append(result)


if __name__ == '__main__':
    # sampling
    tau_samples = []
    for i in range(n_samples):
        tau_samples.append(min_tau + i * (max_tau - min_tau) / (n_samples - 1))
    print(tau_samples)
    # constants
    CONTACTS_PER_HOST = 8
    MU = 1/5  # recovery probability
    SIGMA = 1/17  # immunity lost probability
    BETA = 0.5  # infection probability
    R = 0.01  # recombination probability per allele
    TAU = 0.001  # mutation probability per allele
    GAMMA = 4  # cross-immunity
    N_LOCI = 2
    N_NODES = 300
    SEEDS_PER_STRAIN = 1
    RANDOMNESS = 1  # host contact network randomness, edge reconnecting probability
    N_STEPS = 2000

    # multiprocessing
    result_list = []
    pool = multiprocessing.Pool(processes=4)
    for tau in tau_samples:
        pool.apply_async(simulate, args=(CONTACTS_PER_HOST, MU, SIGMA, BETA, R, tau, GAMMA, N_LOCI, N_NODES, SEEDS_PER_STRAIN, RANDOMNESS, N_STEPS), callback=lambda result: log_result(result, result_list))
    pool.close()
    pool.join()

    # plot
    result_list.sort(key=lambda d: d['tau'])
    plt.figure(figsize=(4, 2))
    diversity = [d['mean_diversity'] for d in result_list]
    discordance = [d['mean_discordance'] for d in result_list]
    dur_list = [d['tau'] for d in result_list]
    plt.plot(dur_list, diversity, 'o-', label='diversity')
    plt.plot(dur_list, discordance, 's--', label='discordance')
    plt.xlim(min_tau, max_tau)
    plt.ylim(0, 1)
    plt.legend()
    plt.xlabel('mutation probability')
    plt.show()
