import networkx as nx
from random import sample, random
import matplotlib.pyplot as plt

n_steps = 50
n_nodes = 500
k_nearest = 4
p_reconnect = 1
beta = 0.3  # infection probability
mu = 0.3  # recovery probability
n_seeds = 4


def check_recovery():
    if random() < mu:
        return 'r'
    else:
        return 'i'


def check_infection():
    if random() < beta:
        return 'i'
    else:
        return 's'


if __name__ == '__main__':
    # generate host contact network
    host_nw = nx.connected_watts_strogatz_graph(n_nodes, k_nearest, p_reconnect)
    # initial setting
    nx.set_node_attributes(host_nw, 's', 'status')
    nx.set_node_attributes(host_nw, 's', 'next_status')
    # seeding
    for n in sample(host_nw.nodes, n_seeds):
        host_nw.nodes[n]['status'] = 'i'
        host_nw.nodes[n]['next_status'] = 'i'
    print(host_nw.nodes.data())
    S = [n_nodes - n_seeds]
    I = [n_seeds]
    R = [0]
    for t in range(1, n_steps):
        # get the infectious nodes
        infectious_nodes = [n for n, v in host_nw.nodes.items() if v['status'] == 'i']
        # infectious nodes recover or/and infect others
        for i in infectious_nodes:
            host_nw.nodes[i]['next_status'] = check_recovery()
            # get the neighbour susceptible nodes
            susceptible_nodes = [n for n in host_nw.adj[i] if host_nw.nodes[n]['status'] == 's']
            # print(susceptible_nodes)
            for s in susceptible_nodes:
                host_nw.nodes[s]['next_status'] = check_infection()
        # apply the status update simultaneously
        for n, v in host_nw.nodes.items():
            v['status'] = v['next_status']
        # record S, I, R
        susceptible_nodes = [n for n, v in host_nw.nodes.items() if v['status'] == 's']
        infectious_nodes = [n for n, v in host_nw.nodes.items() if v['status'] == 'i']
        recovered_nodes = [n for n, v in host_nw.nodes.items() if v['status'] == 'r']
        S.append(len(susceptible_nodes))
        I.append(len(infectious_nodes))
        R.append(len(recovered_nodes))
    # plot
    steps = range(n_steps)
    plt.plot(steps, S, 'b', label='S')
    plt.plot(steps, I, 'r', label='I')
    plt.plot(steps, R, 'g', label='R')
    plt.legend()
    plt.xlabel('Time steps')
    plt.ylabel('Number of hosts')
    plt.title('SIR model')
    plt.show()
