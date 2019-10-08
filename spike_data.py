import json
with open('results_random.txt') as infile:
    results_random = json.load(infile)
results_random = results_random['ran_data']
new_results_random = []
for d in results_random:
    new_d = {'contacts_per_host': d['CONTACTS_PER_HOST'], 'sigma': d['SIGMA'], 'beta': d['BETA'], 'r': d['R'],
             'tao': d['TAO'], 'gamma': d['GAMMA'], 'n_loci': d['N_LOCI'], 'n_nodes': d['N_NODES'], 'mu': d['MU'],
             'mean_diversity': d['RAN_DIV'], 'mean_discordance': d['RAN_DISC']}
    new_results_random.append(new_d)
with open('spike_random.json', 'w') as outfile:
    json.dump(new_results_random, outfile, indent=2)

with open('results_re.txt') as infile:
    results_re = json.load(infile)
results_re = results_re['re_data']
new_results_re = []
for d in results_re:
    new_d = {'contacts_per_host': d['CONTACTS_PER_HOST'], 'sigma': d['SIGMA'], 'beta': d['BETA'], 'r': d['R'],
             'tao': d['TAO'], 'gamma': d['GAMMA'], 'n_loci': d['N_LOCI'], 'n_nodes': d['N_NODES'], 'mu': d['MU'],
             'mean_diversity': d['RE_DIV'], 'mean_discordance': d['RE_DISC']}
    new_results_re.append(new_d)
with open('spike_re.json', 'w') as outfile:
    json.dump(new_results_re, outfile, indent=2)
