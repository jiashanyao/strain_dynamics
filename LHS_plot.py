import matplotlib.pyplot as plt
import json
from random import shuffle

# define scatter color level and inspected parameter
parameter = 'contacts_per_host'   # mu, r

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
result = []
for i in range(len(ran_result)):
    result.append({'ran': ran_result[i], 're':re_result[i]})
shuffle(result)  # shuffle the list to minimize the impact of scatter overlay problem
parameter_list = ['contacts_per_host', 'mu', 'sigma', 'beta', 'r', 'tao', 'gamma', 'n_loci', 'n_nodes']
for parameter in parameter_list:
    plt.figure(figsize=(8, 4))
    plt.subplot(121)
    ran_div = [d['ran']['mean_diversity'] for d in result]
    re_div = [d['re']['mean_diversity'] for d in result]
    colormap = [d['ran'][parameter] for d in result]
    plt.scatter(ran_div, re_div, c=colormap, marker='o', s=1, cmap='viridis')
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
    plt.scatter(ran_disc, re_disc, c=colormap, marker='o', s=1, cmap='viridis')
    plt.plot([0.4, 1], [0.4, 1])
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
    plt.savefig(path + 'LHS_mu=0.1-0.5_r=0.001-0.1_gamma=0.02-10' + parameter + '.png')
    plt.savefig(path + 'LHS_mu=0.1-0.5_r=0.001-0.1_gamma=0.02-10' + parameter + '.pdf')
