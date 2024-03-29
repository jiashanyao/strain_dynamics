import matplotlib.pyplot as plt
import json

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
result = []
for i in range(len(ran_result)):
    result.append({'ran': ran_result[i], 're': re_result[i]})
print(len(result))

plt.figure(figsize=(6, 3))
plt.subplot(121)
ran_div = [d['ran']['mean_diversity'] for d in result]
re_div = [d['re']['mean_diversity'] for d in result]
plt.scatter(ran_div, re_div, marker='o', s=1)
plt.plot([0, 1], [0, 1])
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.xlabel('Diversity in random model')
plt.ylabel('Diversity in regular model')
plt.gca().set_aspect('equal', adjustable='box')
plt.subplot(122)
ran_disc = [d['ran']['mean_discordance'] for d in result]
re_disc = [d['re']['mean_discordance'] for d in result]
plt.scatter(ran_disc, re_disc, marker='o', s=1)
plt.plot([0, 1], [0, 1])
plt.xlim(0.4, 1)
plt.ylim(0.4, 1)
plt.xlabel('Discordance in random model')
plt.ylabel('Discordance in regular model')
plt.gca().set_aspect('equal', adjustable='box')
plt.tight_layout()
# plt.show()
path = 'C:/Users/jiashany/Dropbox/unimelb/ComputingProject/report/images/'
plt.show()
# plt.savefig(path + 'LHS_gamma0.02-5_r0-0.1' + parameter + '.png')
# plt.savefig(path + 'LHS_gamma0.02-5_r0-0.1' + parameter + '.pdf')
