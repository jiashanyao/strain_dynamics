import networkx as nx
from random import sample, random, choice
import matplotlib.pyplot as plt
import math

# constants
N_STEPS = 2000
N_NODES = 250
K_NEAREST = 12
P_RECONNECT = 1
BETA = 0.25  # infection probability
GAMMA = 2  # cross-immunity
MU = 0.14  # recovery probability
SIGMA = 0.043  # immunity lost probability
R = 0.1     # recombination probability per allele
TAO = 0.1     # mutation probability per allele
N_LOCI = 2
N_SEEDS = 16


def generate_strain_space(n_loci):
    strain_space = ['']
    for i in range(n_loci):
        new_strain_space = []
        for s in strain_space:
            new_strain_space.append(s + '0')
            new_strain_space.append(s + '1')
        strain_space = new_strain_space
    return strain_space


def check_immunity_lost(node, current_memory):
    lost = set()
    for strain in current_memory:
        if random() < SIGMA:
            lost.add(strain)
            # print('node', node, 'loses', strain)
    return lost


def check_recovery(node, current_infection):
    recovered = set()
    for strain in current_infection:
        if random() < MU:
            recovered.add(strain)
            # print('node', node, 'recovered from', strain)
    return recovered


def calc_bit_fraction(memory, infect_strain):
    i = 0
    identical_bits = 0
    for bit in infect_strain:
        for strain in memory:
            if bit == strain[i]:
                identical_bits += 1
                break
        i += 1
    return identical_bits / N_LOCI


def recombine(strain1, strain2):
    char_list1 = list(strain1)
    char_list2 = list(strain2)
    for i in range(N_LOCI):
        if random() < R:
            tmp = char_list1[i]
            char_list1[i] = char_list2[i]
            char_list2[i] = tmp
    strain1 = ''.join(char_list1)
    strain2 = ''.join(char_list2)
    return strain1, strain2


def mutate(allele):
    if allele == '0':
        return '1'
    else:
        return '0'


def check_mutation(strain):
    char_list = list(strain)
    for i in range(N_LOCI):
        if random() < TAO:
            char_list[i] = mutate(char_list[i])
    return ''.join(char_list)


def infect_adjacent(network, node):
    if len(network.nodes[node]['current_infection']) == 0:
        return  # skip node that has no infection
    for adj in network.adj[node]:
        # find the strains that are in node's current_infection but not in adj's current_infection
        infection_diff = network.nodes[node]['current_infection'] - network.nodes[adj]['current_infection']
        # if adj's current_infection has all of strains in node's current_infection, skip this adj
        if len(infection_diff) == 0:
            continue
        # arbitrarily select a strain from the strain difference
        infecting_strain = infection_diff.pop()

        # calculate probability
        adj_memory = network.nodes[adj]['current_memory']
        bit_fraction = calc_bit_fraction(adj_memory, infecting_strain)
        vulnerability = (1 - bit_fraction ** (1 / GAMMA)) ** GAMMA
        p_infect = BETA * vulnerability

        # try infecting
        if random() < p_infect:
            infecting_strain = check_mutation(infecting_strain)
            if len(network.nodes[adj]['current_infection']):
                # if adj has current infection, do possible recombination
                infected_strain = set(network.nodes[adj]['current_infection']).pop()
                recombine_output = recombine(infected_strain, infecting_strain)
                network.nodes[adj]['newly_infected'].update(recombine_output)
            else:
                # else no recombination
                network.nodes[adj]['newly_infected'].add(infecting_strain)
                # print('node', node, 'infects node', adj, 'with', infecting_strain)


def print_network(network):
    for n, v in network.nodes.items():
        print(n, v)


def record_host_immune(network, host_immune):
    for immune_record in host_immune.values():
        immune_record.append(0)
    for n, v in network.nodes.items():
        for strain in v['current_memory']:
            host_immune[strain][-1] += 1
    for immune_record in host_immune.values():
        immune_record[-1] = immune_record[-1] / N_NODES


def calc_diversity(strain_freq):
    entropy = 0
    for freq in strain_freq.values():
        if freq != 0:
            entropy += freq * math.log(1 / freq)
    diversity = entropy / math.log(len(strain_freq))
    return diversity


def hamming_dist(strain1, strain2):
    distance = 0
    for i in range(N_LOCI):
        if strain1[i] != strain2[i]:
            distance += 1
    return distance


def calc_discordance(strain_freq):
    numerator = 0
    denominator = 0
    for strain1, freq1 in strain_freq.items():
        for strain2, freq2 in strain_freq.items():
            if strain1 != strain2:
                numerator += hamming_dist(strain1, strain2) * freq1 * freq2
                denominator += freq1 * freq2
    discordance = numerator / denominator / N_LOCI
    return discordance


def main():
    # generate host contact network
    host_nw = nx.connected_watts_strogatz_graph(N_NODES, K_NEAREST, P_RECONNECT)
    # nx.draw(host_nw, with_labels=True)
    # plt.show()

    # data field setting
    for n, v in host_nw.nodes.items():
        v['current_memory'] = set()     # actual memory
        v['current_infection'] = set()  # actual infected strains
        v['newly_lost'] = set()         # temporary hold for newly lost strain immunity
        v['newly_infected'] = set()     # temporary hold for newly infected strains
        v['newly_recovered'] = set()    # temporary hold for newly recovered strains

    # generate strain space
    strain_space = generate_strain_space(N_LOCI)
    # print('strain space:', strain_space)

    # initialize population record
    host_immune = {}
    for strain in strain_space:
        host_immune[strain] = []

    # seeding
    for n in sample(host_nw.nodes, N_SEEDS):
        seeded_strain = choice(strain_space)
        host_nw.nodes[n]['current_infection'].add(seeded_strain)
        print('seeded', seeded_strain)

    # print('step 0 (seeding)')
    # print_network(host_nw)

    # record strain population of step 0
    record_host_immune(host_nw, host_immune)

    # simulate
    for t in range(1, N_STEPS):
        # print('step', t)

        # records infection and memory changes in temporary data fields for each node
        for n, v in host_nw.nodes.items():
            v['newly_lost'] = check_immunity_lost(n, v['current_memory'])
            v['newly_recovered'] = check_recovery(n, v['current_infection'])
            infect_adjacent(host_nw, n)

        # print('before update:')
        # print_network(host_nw)

        # apply change records and reset the temporary hold fields
        for n, v in host_nw.nodes.items():
            # apply changes
            v['current_memory'] = (v['current_memory'] | v['newly_recovered']) - v['newly_lost']
            v['current_infection'] = (v['current_infection'] | v['newly_infected']) - v['newly_recovered']
            # reset temporary hold fields
            v['newly_lost'] = set()
            v['newly_recovered'] = set()
            v['newly_infected'] = set()

        # print('after update:')
        # print_network(host_nw)

        record_host_immune(host_nw, host_immune)

    # plot
    steps = range(N_STEPS)
    for strain, host_immune_record in host_immune.items():
        plt.plot(steps, host_immune_record, label=strain)
    plt.legend()
    plt.xlabel('Time steps')
    plt.ylabel('Hosts immune to a strain')
    plt.show()


if __name__ == '__main__':
    main()
