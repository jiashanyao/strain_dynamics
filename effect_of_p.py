import multiprocessing
from matplotlib import pyplot as plt
from multistrain_model import simulate

max_power = -4
n_samples = 10


def log_result(result, result_list):
    if result:
        result_list.append(result)


if __name__ == '__main__':
    # sampling
    p_samples = [1]
    for i in range(n_samples - 1):
        p_samples.append(p_samples[-1] * 10 ** (max_power / (n_samples - 1)))
    print(p_samples)
    # constants
    CONTACTS_PER_HOST = 12
    MU = 1/5  # recovery probability
    SIGMA = 1/20  # immunity lost probability
    BETA = 0.3  # infection probability
    R = 0.01  # recombination probability per allele
    TAU = 0.001  # mutation probability per allele
    GAMMA = 4  # cross-immunity
    N_LOCI = 4
    N_NODES = 500
    SEEDS_PER_STRAIN = 1
    RANDOMNESS = 1  # host contact network randomness, edge reconnecting probability
    N_STEPS = 2000
    parameters = [CONTACTS_PER_HOST, MU, SIGMA, BETA, R, TAU, GAMMA, N_LOCI, N_NODES]

    # multiprocessing
    result_list = []
    pool = multiprocessing.Pool(processes=4)
    for p in p_samples:
        pool.apply_async(simulate, args=(*parameters, SEEDS_PER_STRAIN, p, N_STEPS), callback=lambda result: log_result(result, result_list))
    pool.close()
    pool.join()

    # plot
    result_list.sort(key=lambda d: d['randomness'])
    plt.figure(figsize=(4, 4))
    plt.subplot(211)
    diversity = [d['mean_diversity'] for d in result_list]
    discordance = [d['mean_discordance'] for d in result_list]
    p_list = [d['randomness'] for d in result_list]
    plt.plot(p_list, diversity, 'o-', label='diversity')
    plt.plot(p_list, discordance, 's--', label='discordance')
    plt.xscale('log')
    plt.xlabel('p')
    plt.ylim(0, 1)
    plt.legend()
    plt.subplot(212)
    cluster_coefficient = [d['cluster_coefficient'] for d in result_list]
    plt.plot(p_list, cluster_coefficient, 'o-', label='clustering coefficient')
    plt.xscale('log')
    plt.xlabel('p')
    plt.ylim(0, 1)
    plt.legend()
    plt.show()
