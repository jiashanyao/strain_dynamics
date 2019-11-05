import multiprocessing
from matplotlib import pyplot as plt
from multistrain_model import simulate

min_beta = 0.2
max_beta = 0.8
n_samples = 13


def log_result(result, result_list):
    if result:
        result_list.append(result)


if __name__ == '__main__':
    # sampling
    beta_samples = []
    for i in range(n_samples):
        beta_samples.append(min_beta + i * (max_beta - min_beta) / (n_samples - 1))
    print(beta_samples)
    # constants
    CONTACTS_PER_HOST = 8
    MU = 1/6  # recovery probability
    SIGMA = 1/20  # immunity lost probability
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
    for beta in beta_samples:
        pool.apply_async(simulate, args=(CONTACTS_PER_HOST, MU, SIGMA, beta, R, TAU, GAMMA, N_LOCI, N_NODES, SEEDS_PER_STRAIN, RANDOMNESS, N_STEPS), callback=lambda result: log_result(result, result_list))
    pool.close()
    pool.join()

    # plot
    result_list.sort(key=lambda d: d['beta'])
    plt.figure(figsize=(4, 2))
    diversity = [d['mean_diversity'] for d in result_list]
    discordance = [d['mean_discordance'] for d in result_list]
    dur_list = [d['beta'] for d in result_list]
    plt.plot(dur_list, diversity, 'o-', label='diversity')
    plt.plot(dur_list, discordance, 's--', label='discordance')
    plt.xlim(min_beta, max_beta)
    plt.ylim(0, 1)
    plt.legend()
    plt.xlabel('Infection probability')
    plt.tight_layout()
    plt.show()
