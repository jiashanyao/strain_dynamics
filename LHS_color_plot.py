import matplotlib.pyplot as plt
import json
from random import shuffle

with open('ran_new_recombination.json') as infile:
    ran_result = json.load(infile)
with open('re_new_recombination.json') as infile:
    re_result = json.load(infile)

# pre-process the data: trim uncommon entries and sort the trim results to match ran_result and re_result
ran_beta_keys = set([d['beta'] for d in ran_result])
re_beta_keys = set(d['beta'] for d in re_result)
common_beta_keys = ran_beta_keys & re_beta_keys
ran_result = [d for d in ran_result if d['beta'] in common_beta_keys]
re_result = [d for d in re_result if d['beta'] in common_beta_keys]
ran_result.sort(key=lambda d: d['beta'])
re_result.sort(key=lambda d: d['beta'])
result = []
for i in range(len(ran_result)):
    result.append({'ran': ran_result[i], 're': re_result[i]})
shuffle(result)  # shuffle the list to minimize the impact of scatter overlay problem
parameter_list = ['contacts_per_host', 'mu', 'sigma', 'beta', 'r', 'gamma', 'n_loci', 'n_nodes']
for parameter in parameter_list:

    plt.figure(figsize=(8, 4))
    plt.subplot(121)
    ran_div = [d['ran']['mean_diversity'] for d in result]
    re_div = [d['re']['mean_diversity'] for d in result]
    colormap = [d['ran'][parameter] for d in result]
    if parameter in ['mu', 'sigma']:
        for i in range(len(colormap)):
            colormap[i] = 1/colormap[i]
    plt.scatter(ran_div, re_div, c=colormap, marker='o', s=1)
    plt.plot([0, 1], [0, 1])
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.colorbar()
    plt.xlabel('Diversity in random model')
    plt.ylabel('Diversity in regular model')
    plt.title(parameter + ' in different color')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.subplot(122)
    ran_disc = [d['ran']['mean_discordance'] for d in result]
    re_disc = [d['re']['mean_discordance'] for d in result]
    colormap = [d['ran'][parameter] for d in result]
    if parameter in ['mu', 'sigma']:
        for i in range(len(colormap)):
            colormap[i] = 1/colormap[i]
    plt.scatter(ran_disc, re_disc, c=colormap, marker='o', s=1)
    plt.plot([0, 1], [0, 1])
    plt.xlim(0.4, 1)
    plt.ylim(0.4, 1)
    plt.colorbar()
    plt.xlabel('Discordance in random model')
    plt.ylabel('Discordance in regular model')
    plt.title(parameter + ' in different color')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.tight_layout()
    # plt.show()
    path = 'C:/Users/jiashany/Dropbox/unimelb/ComputingProject/report/images/'
    plt.show()
    # plt.savefig(path + 'LHS_gamma0.02-5_r0-0.1' + parameter + '.png')
    # plt.savefig(path + 'LHS_gamma0.02-5_r0-0.1' + parameter + '.pdf')
