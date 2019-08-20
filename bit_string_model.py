import networkx as nx
from random import sample, random, choice
import matplotlib.pyplot as plt

n_steps = 50
n_nodes = 10
k_nearest = 4
p_reconnect = 1
beta = 0.3  # infection probability
gamma = 0.02  # cross-immunity
mu = 0.3  # recovery probability
sigma = 0.1  # immunity lost probability
n_loci = 3
n_seeds = 4


def generate_strain_space(n_loci):
    strain_space = ['']
    for i in range(n_loci):
        new_strain_space = []
        for s in strain_space:
            new_strain_space.append(s + '0')
            new_strain_space.append(s + '1')
        strain_space = new_strain_space
    return strain_space


def check_immunity_lost(current_memory):
    lost = set()
    for strain in current_memory:
        if random() < sigma:
            lost.add(strain)
    return lost


def check_recovery(current_infection):
    recovered = set()
    for strain in current_infection:
        if random() < mu:
            recovered.add(strain)
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
    return identical_bits/n_loci


def infect_adjacent(network, node):
    for adj in network.adj[node]:
        # find the strains that are in node's current_infection but not in adj's current_infection
        infection_diff = network.nodes[node]['current_infection'] - network.nodes[adj]['current_infection']

        # if adj's current_infection has all of strains in node's current_infection, skip this adj
        if len(infection_diff) == 0:
            continue

        # arbitrarily select a strain from the strain difference
        infect_strain = infection_diff.pop()

        # calculate probability
        adj_memory = network.nodes[adj]['current_memory']
        bit_fraction = calc_bit_fraction(adj_memory, infect_strain)
        vulnerability = (1 - bit_fraction ** (1/gamma)) ** gamma
        p_infect = beta * vulnerability

        # start infecting
        if random() < p_infect:
            # add infecting strain to adj's newly_infected set
            network.nodes[adj]['newly_infected'].add(infect_strain)


if __name__ == '__main__':
    # generate host contact network
    host_nw = nx.connected_watts_strogatz_graph(n_nodes, k_nearest, p_reconnect)

    # data structure setting
    for n, v in host_nw.nodes.items():

        # data fields for actual host properties
        v['current_memory'] = set()
        v['current_infection'] = set()

        # data fields for recording changes in host properties
        v['next_memory'] = set()
        v['newly_infected'] = set()
        v['newly_recovered'] = set()

    # generate strain space
    strain_space = generate_strain_space(n_loci)
    # print(strain_space)

    # seeding
    for n in sample(host_nw.nodes, n_seeds):
        seeded_strain = choice(strain_space)
        host_nw.nodes[n]['current_infection'].add(seeded_strain)
        host_nw.nodes[n]['current_memory'].add(seeded_strain)
    print(host_nw.nodes.data())

    # simulate
    for t in range(n_steps):

        # records infection and memory changes in temporary data fields for each node
        for n, v in host_nw.nodes.items():

            # skip uninfected nodes
            if not v['current_infection']:
                continue

            # check immunity lost
            v['next_memory'] = check_immunity_lost(v['current_memory'])

            # check recovery
            v['newly_recovered'] = check_recovery(v['current_infection'])

            # check infection
            infect_adjacent(host_nw, n)

        # apply change records and reset the temporary hold fields
        for n, v in host_nw.nodes.items():

            # apply changes
            v['current_memory'] = v['next_memory']  # shallow copy as next_memory will be assigned with a new object
            v['current_infection'] = v['current_infection'] - v['newly_recovered']
            v['current_infection'] = v['current_infection'] | v['newly_infected']

            # reset temporary hold fields
            v['new_memory'] = set()
            v['newly_recovered'] = set()
            v['newly_infected'] = set()
