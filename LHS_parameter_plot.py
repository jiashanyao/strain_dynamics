import json
import matplotlib.pyplot as plt

with open('ran_gamma=0.02-20_u=2-10_c=8-12.json') as infile:
    ran_result = json.load(infile)
with open('re_gamma=0.02-20_u=2-10_c=8-12.json') as infile:
    re_result = json.load(infile)

# pre-process the data: trim uncommon entries and sort the trim results to match ran_result and re_result
ran_beta_keys = set([d['beta'] for d in ran_result])
re_beta_keys = set(d['beta'] for d in re_result)
common_beta_keys = ran_beta_keys & re_beta_keys
ran_result = [d for d in ran_result if d['beta'] in common_beta_keys]
re_result = [d for d in re_result if d['beta'] in common_beta_keys]
ran_result.sort(key=lambda d: d['beta'])
re_result.sort(key=lambda d: d['beta'])

parameter_list = ['mu', 'sigma']
path = 'C:/Users/jiashany/Dropbox/unimelb/ComputingProject/report/images/'
colormap = [d['n_loci'] for d in ran_result]
for parameter in parameter_list:
    plt.figure(figsize=(7, 3))
    plt.subplot(121)
    parameter_value = [1/d[parameter] for d in ran_result]
    diversity = [d['mean_diversity'] for d in ran_result]
    plt.scatter(parameter_value, diversity, marker='o', s=1)
    plt.ylim(0, 1)
    plt.xlabel('1/' + parameter)
    plt.ylabel('diversity random')
    plt.subplot(122)
    parameter_value = [1/d[parameter] for d in re_result]
    diversity = [d['mean_diversity'] for d in re_result]
    plt.scatter(parameter_value, diversity, marker='o', s=1)
    plt.ylim(0, 1)
    plt.xlabel('1/' + parameter)
    plt.ylabel('diversity regular')
    plt.tight_layout()
    # plt.savefig(path + 'diversity_' + parameter + '.pdf')
    # plt.savefig(path + 'diversity_' + parameter + '.png')
    plt.show()

    plt.figure(figsize=(7, 3))
    plt.subplot(121)
    parameter_value = [1/d[parameter] for d in ran_result]
    discordance = [d['mean_discordance'] for d in ran_result]
    plt.scatter(parameter_value, discordance, marker='o', s=1)
    plt.ylim(0, 1)
    plt.xlabel('1/' + parameter)
    plt.ylabel('discordance random')
    plt.subplot(122)
    parameter_value = [1/d[parameter] for d in re_result]
    discordance = [d['mean_discordance'] for d in re_result]
    plt.scatter(parameter_value, discordance, marker='o', s=1)
    plt.ylim(0, 1)
    plt.xlabel('1/' + parameter)
    plt.ylabel('discordance regular')
    plt.tight_layout()
    # plt.savefig(path + 'discordance_' + parameter + '.pdf')
    # plt.savefig(path + 'discordance_' + parameter + '.png')
    plt.show()
