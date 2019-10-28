from effect_of_mu_and_sigma import *

with open('result_matrix_beta=0.2.json') as infile:
    result_matrix = json.load(infile)
plt.figure(figsize=(4, 3.5))
plt.imshow(result_matrix, extent=[min_imm_dur, max_imm_dur, max_inf_dur, min_inf_dur], aspect=2.857)
plt.colorbar()
plt.xlabel('Immunity duration')
plt.ylabel('Infection duration')
plt.title('Discordance')
plt.tight_layout()
plt.show()
