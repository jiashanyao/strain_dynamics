import networkx as nx
from random import sample, random, choice
import matplotlib.pyplot as plt
import math


def generate_strain_space(n_loci):
    strain_space = ['']
    for i in range(n_loci):
        new_strain_space = []
        for s in strain_space:
            new_strain_space.append(s + '0')
            new_strain_space.append(s + '1')
        strain_space = new_strain_space
    return strain_space


def check_immunity_lost(node, current_memory, sigma):
    lost = set()
    for strain in current_memory:
        if random() < sigma:
            lost.add(strain)
            # print('node', node, 'loses', strain)
    return lost


def check_recovery(node, current_infection, mu):
    recovered = set()
    for strain in current_infection:
        if random() < mu:
            recovered.add(strain)
            # print('node', node, 'recovered from', strain)
    return recovered


def calc_bit_fraction(memory, infect_strain, n_loci):
    i = 0
    identical_bits = 0
    for bit in infect_strain:
        for strain in memory:
            if bit == strain[i]:
                identical_bits += 1
                break
        i += 1
    return identical_bits / n_loci


def recombine(strain1, strain2, n_loci, r):
    char_list1 = list(strain1)
    char_list2 = list(strain2)
    for i in range(n_loci):
        if random() < r:
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


def check_mutation(strain, n_loci, tao):
    char_list = list(strain)
    for i in range(n_loci):
        if random() < tao:
            char_list[i] = mutate(char_list[i])
    return ''.join(char_list)


def infect_adjacent(network, node, n_loci, gamma, beta, tao, r):
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
        bit_fraction = calc_bit_fraction(adj_memory, infecting_strain, n_loci)
        vulnerability = (1 - bit_fraction ** (1 / gamma)) ** gamma
        p_infect = beta * vulnerability

        # try infecting
        if random() < p_infect:
            infecting_strain = check_mutation(infecting_strain, n_loci, tao)
            if len(network.nodes[adj]['current_infection']):
                # if adj has current infection, do possible recombination
                infected_strain = set(network.nodes[adj]['current_infection']).pop()
                recombine_output = recombine(infected_strain, infecting_strain, n_loci, r)
                network.nodes[adj]['newly_infected'].update(recombine_output)
            else:
                # else no recombination
                network.nodes[adj]['newly_infected'].add(infecting_strain)
                # print('node', node, 'infects node', adj, 'with', infecting_strain)


def print_network(network):
    for n, v in network.nodes.items():
        print(n, v)


def record_host_immune(network, host_immune, n_nodes):
    for immune_record in host_immune.values():
        immune_record.append(0)
    for n, v in network.nodes.items():
        for strain in v['current_memory']:
            host_immune[strain][-1] += 1
    for immune_record in host_immune.values():
        immune_record[-1] = immune_record[-1] / n_nodes


def calc_diversity(strain_freq):
    entropy = 0
    for freq in strain_freq.values():
        if freq != 0:
            entropy += freq * math.log(1 / freq)
    diversity = entropy / math.log(len(strain_freq))
    return diversity


def hamming_dist(strain1, strain2, n_loci):
    distance = 0
    for i in range(n_loci):
        if strain1[i] != strain2[i]:
            distance += 1
    return distance


def calc_discordance(strain_freq, n_loci):
    numerator = 0
    denominator = 0
    for strain1, freq1 in strain_freq.items():
        for strain2, freq2 in strain_freq.items():
            if strain1 != strain2:
                numerator += hamming_dist(strain1, strain2, n_loci) * freq1 * freq2
                denominator += freq1 * freq2
    discordance = numerator / denominator / n_loci
    return discordance


def seed(network, n_seeds, strain):
    candidates = [n for n, v in network.nodes.items() if strain not in v['current_infection']]
    sampled_nodes = sample(candidates, n_seeds)
    for node in sampled_nodes:
        network.nodes[node]['current_infection'].add(strain)
        print('seeded', node, 'with', strain)


def simulate(contacts_per_host, mu, sigma, beta, r, tao, gamma, n_loci, n_nodes, randomness, n_steps, seed_sequence=None):
    print('contacts_per_host', contacts_per_host)
    print('mu', mu)
    print('sigma', sigma)
    print('beta', beta)
    print('r', r)
    print('tao', tao)
    print('gamma', gamma)
    print('n_loci', n_loci)
    print('n_nodes', n_nodes)
    print('randomness', randomness)
    print('n_steps', n_steps)

    # generate host contact network
    host_nw = nx.connected_watts_strogatz_graph(n_nodes, contacts_per_host, randomness)

    # data field setting
    for n, v in host_nw.nodes.items():
        v['current_memory'] = set()     # actual memory
        v['current_infection'] = set()  # actual infected strains
        v['newly_lost'] = set()         # temporary hold for newly lost strain immunity
        v['newly_infected'] = set()     # temporary hold for newly infected strains
        v['newly_recovered'] = set()    # temporary hold for newly recovered strains

    # if seed_sequence is provided, overwrite the n_loci with the provided strain length
    if seed_sequence:
        n_loci = len(seed_sequence[0][1])

    # generate strain space
    strain_space = generate_strain_space(n_loci)

    # initialize population record
    host_immune = {}
    for strain in strain_space:
        host_immune[strain] = []

    # default seeding if seed_sequence is not provided
    if not seed_sequence:
        # seed each strain in strain space to a host
        for strain_seed in strain_space:
            seed(host_nw, 1, strain_seed)

    # print('step 0')
    # print_network(host_nw)

    # record strain population of step 0
    record_host_immune(host_nw, host_immune, n_nodes)

    # simulate
    for t in range(1, n_steps):
        # print('step', t)

        # if in manual mode, seed another strain to the community at some time step
        if seed_sequence:
            if t == seed_sequence[0][0]:
                seed(host_nw, 1, seed_sequence[0][1])
                seed_sequence.pop(0)

        # records infection and memory changes in temporary data fields for each node
        for n, v in host_nw.nodes.items():
            v['newly_lost'] = check_immunity_lost(n, v['current_memory'], sigma)
            v['newly_recovered'] = check_recovery(n, v['current_infection'], mu)
            infect_adjacent(host_nw, n, n_loci, gamma, beta, tao, r)

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

        record_host_immune(host_nw, host_immune, n_nodes)

    # plot
    steps = range(n_steps)
    for strain, host_immune_record in host_immune.items():
        plt.plot(steps, host_immune_record, label=strain)
    plt.legend()
    plt.xlabel('Time steps')
    plt.ylabel('Hosts immune to a strain')
    plt.show()


if __name__ == '__main__':
    # constants
    CONTACTS_PER_HOST = 10
    MU = 0.14  # recovery probability
    SIGMA = 0.03  # immunity lost probability
    BETA = 0.2  # infection probability
    R = 0.0  # recombination probability per allele
    TAO = 0.000  # mutation probability per allele
    GAMMA = 3  # cross-immunity
    N_LOCI = 3
    N_NODES = 200
    RANDOMNESS = 1  # host contact network randomness, edge reconnecting probability
    N_STEPS = 2000
    parameters = [CONTACTS_PER_HOST, MU, SIGMA, BETA, R, TAO, GAMMA, N_LOCI, N_NODES]
    simulate(*parameters, RANDOMNESS, N_STEPS, [[1,'00'], [200,'01']])
