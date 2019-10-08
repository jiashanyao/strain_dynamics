import multiprocessing

import matplotlib.pyplot as plt
import json
from multistrain_model import simulate

if __name__ == '__main__':

    with open('ran_mu=0.1-0.5_r=0.001-0.1_gamma=0.02-10.json') as infile:
        ran_result = json.load(infile)
    with open('re_mu=0.1-0.5_r=0.001-0.1_gamma=0.02-10.json') as infile:
        re_result = json.load(infile)

    # pre-process the data: trim uncommon entries and sort the trim results to match ran_result and re_result
    ran_beta_keys = set([d['beta'] for d in ran_result])
    re_beta_keys = set(d['beta'] for d in re_result)
    common_beta_keys = ran_beta_keys & re_beta_keys
    ran_result = [d for d in ran_result if d['beta'] in common_beta_keys]
    re_result = [d for d in re_result if d['beta'] in common_beta_keys]
    ran_result.sort(key=lambda d: d['beta'])
    re_result.sort(key=lambda d: d['beta'])

    pool = multiprocessing.Pool()

    for i in range(len(ran_result)):
        if 0.65 < ran_result[i]['mean_discordance'] < 0.7 and 0.65 < re_result[i]['mean_discordance'] < 0.7:
            result = ran_result[i]
            pool.apply_async(simulate, args=(result['contacts_per_host'], result['mu'], result['sigma'], result['beta'], result['r'], result['tao'],
                                             result['gamma'], result['n_loci'], result['n_nodes'], result['seeds_per_strain'], result['randomness'],
                                             result['n_steps'], None, True, True))
    pool.close()
    pool.join()
