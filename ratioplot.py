from effect_of_mu_and_sigma import *
from matplotlib import pyplot as plt

with open('result_matrix_beta=0.2.json') as infile:
    result_matrix = json.load(infile)

ratio_list = []
discordance_list = []
for row in range(n):
    for col in range(n):
        inf_dur = min_inf_dur + row/n * (max_inf_dur - min_inf_dur)
        imm_dur = min_imm_dur + col/n * (max_imm_dur - min_imm_dur)
        ratio_list.append(imm_dur/inf_dur)
        discordance_list.append(result_matrix[row][col])
plt.figure(figsize=(4, 2))
plt.scatter(ratio_list, discordance_list, marker='o', s=1)
plt.xlabel('Immunity duration / Infection duration')
plt.ylabel('Discordance')
plt.ylim(0.65, 1)
plt.tight_layout()
plt.show()
