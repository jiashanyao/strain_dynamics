import multiprocessing
from multistrain_model import simulate
import json

min_inf_dur = 3
max_inf_dur = 10
min_imm_dur = 10
max_imm_dur = 30
n = 50


def log_result(result, result_matrix):
    if result:
        row = result['row']
        col = result['col']
        print('----------------------', row, col, '-------------------------')
        result_matrix[row][col] = result['mean_discordance']


if __name__ == '__main__':
    # constants
    CONTACTS_PER_HOST = 8
    MU = 1/6  # recovery probability
    SIGMA = 1/20  # immunity lost probability
    BETA = 0.2  # infection probability
    R = 0.01  # recombination probability per allele
    TAU = 0.001  # mutation probability per allele
    GAMMA = 4  # cross-immunity
    N_LOCI = 2
    N_NODES = 300
    SEEDS_PER_STRAIN = 1
    RANDOMNESS = 1  # host contact network randomness, edge reconnecting probability
    N_STEPS = 2000

    # multiprocessing
    result_matrix = [[0 for col in range(n)] for row in range(n)]
    pool = multiprocessing.Pool(processes=4)
    for i in range(n):
        inf_dur = min_inf_dur + i * (max_inf_dur - min_inf_dur) / (n - 1)
        for j in range(n):
            imm_dur = min_imm_dur + j * (max_imm_dur - min_imm_dur) / (n - 1)
            pool.apply_async(simulate, args=(CONTACTS_PER_HOST, 1/inf_dur, 1/imm_dur, BETA, R, TAU, GAMMA,
                                             N_LOCI, N_NODES, SEEDS_PER_STRAIN, RANDOMNESS, N_STEPS,
                                             None, False, False, i, j),
                             callback=lambda result: log_result(result, result_matrix))
    pool.close()
    pool.join()
    with open('result_matrix_beta=0.2.json', 'w') as outfile:
        json.dump(result_matrix, outfile, indent=2)
