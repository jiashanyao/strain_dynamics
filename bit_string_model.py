import networkx as nx
from random import sample, random, choice
import matplotlib.pyplot as plt

# constants
N_STEPS = 500
N_NODES = 500
K_NEAREST = 12
P_RECONNECT = 0
BETA = 0.25  # infection probability
GAMMA = 2  # cross-immunity
MU = 0.14  # recovery probability
SIGMA = 0.043  # immunity lost probability
N_LOCI = 2
N_SEEDS = 50


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
        vulnerability = (1 - bit_fraction ** (1 / GAMMA)) ** GAMMA
        p_infect = BETA * vulnerability

        # start infecting
        if random() < p_infect:
            # add infecting strain to adj's newly_infected set
            network.nodes[adj]['newly_infected'].add(infect_strain)
            # print('node', node, 'infects node', adj, 'with', infect_strain)


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


def main():
    # generate host contact network
    host_nw = nx.connected_watts_strogatz_graph(N_NODES, K_NEAREST, P_RECONNECT)
    # nx.draw(host_nw, with_labels=True)
    # plt.show()

    # data structure setting
    for n, v in host_nw.nodes.items():

        # data fields for actual host properties
        v['current_memory'] = set()
        v['current_infection'] = set()

        # data fields for recording changes in host properties
        v['newly_lost'] = set()
        v['newly_infected'] = set()
        v['newly_recovered'] = set()

    # generate strain space
    strain_space = generate_strain_space(N_LOCI)
    print('strain space:', strain_space)

    # initialize population record
    host_immune = {}
    for strain in strain_space:
        host_immune[strain] = []

    # seeding
    for n in sample(host_nw.nodes, N_SEEDS):
        seeded_strain = choice(strain_space)
        host_nw.nodes[n]['current_infection'].add(seeded_strain)
    # print('step 0 (seeding)')
    # print_network(host_nw)

    # record strain population of step 0
    record_host_immune(host_nw, host_immune)
    print(host_immune)

    # simulate
    for t in range(1, N_STEPS):
        # print('step', t)

        # records infection and memory changes in temporary data fields for each node
        for n, v in host_nw.nodes.items():

            # skip uninfected nodes
            if not v['current_infection']:
                continue

            # check immunity lost
            v['newly_lost'] = check_immunity_lost(n, v['current_memory'])

            # check recovery
            v['newly_recovered'] = check_recovery(n, v['current_infection'])

            # check infection
            infect_adjacent(host_nw, n)

        # print('before update:')
        # print_network(host_nw)

        # apply change records and reset the temporary hold fields
        for n, v in host_nw.nodes.items():

            # apply changes
            v['current_memory'] = (v['current_memory'] - v['newly_lost']) | v['newly_recovered']
            v['current_infection'] = (v['current_infection'] - v['newly_recovered']) | v['newly_infected']

            # reset temporary hold fields
            v['newly_lost'] = set()
            v['newly_recovered'] = set()
            v['newly_infected'] = set()

        # print('after update:')
        # print_network(host_nw)

        record_host_immune(host_nw, host_immune)
        # print(population)

    # calculate
    steps = range(N_STEPS)
    for strain, host_immune_record in host_immune.items():
        plt.plot(steps, host_immune_record, label=strain)
    plt.legend()
    plt.xlabel('Time steps')
    plt.ylabel('Hosts immune to a strain')
    plt.show()


if __name__ == '__main__':
    main()
