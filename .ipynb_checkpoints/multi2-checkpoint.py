from multiprocessing import Pool, Manager, Process
import networkx as nx
from networkx.algorithms.coloring.greedy_coloring import greedy_color
import numpy as np
import itertools
import time

def node_chunks(G, p):
    ''' Given a graph G of n nodes, and an integer p, split the set of vertices into p sets of 
    length |V_i| such that floor(n/p) <= |V_i| <= ceiling(n/p) '''
    k, m = divmod(len(G.nodes), p)
    chunks = [list(G.nodes)[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(p)]
    return chunks


def induce_subgraph(G, u, block):
    ''' Given a block of vertices that will be colored in parallel, and a vertex u in that
    block, construct a graph induced on the vertex set V(G) - {block} + {u} '''
    nodes = set(G.nodes) - set(block)
    nodes.add(u)
    return G.subgraph(nodes)


def nodes_iter(Vp):
    ''' Given p chunks of nodes, order the vertices into an array, blocks, such that blocks[i]
    indicates vertices that will be colored in parallel at iteration i '''
    batches = list(np.array(Vp).T)
    return batches


def color_node(G, u, coloring):
    ''' Given a vertex, color it based on the graph and previously assigned coloring to other
    vertices '''
    neighbour_colors = {coloring[v] for v in G[u] if v in coloring}
    for color in itertools.count():
        if color not in neighbour_colors:
            break
    coloring[u] = color
    return (u, color)


def parallel_color(G, p):
    ''' Implementation of algorithm 1 from the research paper with a slight modification:
        The Pool class in Python uses a method of interleaving for function calls, so
        true parallelism cannot be acheived in the same way of the paper.
        
        Instead, at each iteration in parallel, color_node is called with the vertex set
        of G.nodes - {other vertices to be colored in this parallel iteration}, and the
        induced subgraph on this set. Thus, colors that are assigned that would still be
        seen due to interleaving are not yet assigned. '''
    # Phase 1
    Vp = node_chunks(G, p)
    min_blen = min([len(block) for block in Vp])
    batches = nodes_iter(Vp)

    manager = Manager()
    coloring = manager.dict()

    for i in range(len(batches)):

        jobs = [Process(target=color_node, args=(G, batches[i][j], coloring,)) for j in range(p)]
        _ = [proc.start() for proc in jobs]
        _ = [proc.join() for proc in jobs]
    
    # Phase 2
    A = []

    for i in range(len(batches)):
        for v in batches[i]:
            for u in G.neighbors(v):
                if u in batches[i]:
                    if coloring[u] == coloring[v]:
                        if min(u, v) not in A:
                            A.append(min(u, v))

    for v in A:
        color_node(G, v, coloring)
    return coloring, len(A)


def parallel_color_strat(G, p, blocking):
    ''' Modified Algorithm 1 to take alternate blocking algorithms.'''
    # Phase 1
    # Takes blocking function 
    Vp = blocking(G, p)
    min_blen = min([len(block) for block in Vp])
    batches = nodes_iter(Vp)

    manager = Manager()
    coloring = manager.dict()

    for i in range(len(batches)):
        batches[i][1]

        jobs = [Process(target=color_node, args=(G, batches[i][j], coloring,)) for j in range(p)]
        _ = [proc.start() for proc in jobs]
        _ = [proc.join() for proc in jobs]
    
    # Phase 2
    A = []

    for i in range(len(batches)):
        for v in batches[i]:
            for u in G.neighbors(v):
                if u in batches[i]:
                    if coloring[u] == coloring[v]:
                        if min(u, v) not in A:
                            A.append(min(u, v))

    for v in A:
        color_node(G, v, coloring)
    return coloring, len(A)


def avg_degree(G):
    ''' Find the average degree of a graph, G '''
    sum_d = 0
    for u in G.nodes:
        sum_d += G.degree(u)
    return sum_d / len(G.nodes)


def random(G, p):
    k, m = divmod(len(G.nodes), p)
    chunks = [list(G.nodes)[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(p)]
    return chunks


def descending_degree_chunks(G, p):
    sorted_nodes = sorted(G, key=G.degree, reverse=True)
    k, m = divmod(len(sorted_nodes), p)
    chunks = [sorted_nodes[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(p)]
    return chunks


def descending_degree_within(G, p):
    sorted_nodes = sorted(G, key=G.degree, reverse=True)
    print(sorted_nodes)
    chunks = [[] for _ in range(p)]
    while (sorted_nodes):
        for i in range(p):
            if sorted_nodes:
                chunks[i].append(sorted_nodes.pop(0))
    return chunks


if __name__ == "__main__":
    # We will test networkx's default greedy_color algorithm, that colors with a strategy
    # of coloring highest degree vertices first against our parallelized implementation.
    # The respective run-times and chromatic number that the algorithms return will be
    # compared, while varying average degree and number of processors
    # The number of conflicts in our algorithm will be analyzed as well

    # 1 Processor
    n = 100
    prob = .1
    p = 1
    proc1_avg_deg = {}
    proc1_runtime_parallel = {}
    proc1_chrom_num_parallel = {}
    proc1_num_conflicts = {}
    proc1_runtime_greedy = {}
    proc1_chrom_num_greedy = {}


    while n < 1000:
        while prob <= .6:
            G = nx.gnp_random_graph(n, prob)
            #print(
            #    'Using a graph generated from G_', n, ',', prob, 'with average degree', avg_degree(G) 
            #)
            proc1_avg_deg[str(n) + ' ' + str(prob)] = avg_degree(G)

            start_time = time.process_time()
            color_p, num_conflicts = parallel_color(G, p)
            #print(
            #    'The parallelized algorithm took', time.process_time() - start_time,
            #    'seconds to run with', num_conflicts, 'conflicts and a chromatic number of',
            #    max(color_p.values())
            #    )
            proc1_runtime_parallel[str(n) + ' ' + str(prob)] = time.process_time() - start_time
            proc1_num_conflicts[str(n) + ' ' + str(prob)] = num_conflicts
            proc1_chrom_num_parallel[str(n) + ' ' + str(prob)] = max(color_p.values())

            start_time = time.process_time()
            color_g = greedy_color(G)
            #print(
            #    'Networkx greedy algorithm took', time.process_time() - start_time,
            #    'seconds to run with a chromatic number of', max(color_g.values())
            #    )
            proc1_runtime_greedy[str(n) + ' ' + str(prob)] = time.process_time() - start_time
            proc1_chrom_num_greedy = max(color_g.values())

            n += 100
            prob += .1
    
    # 2 Processors
    n = 100
    prob = .1
    p = 2
    proc2_avg_deg = {}
    proc2_runtime_parallel = {}
    proc2_chrom_num_parallel = {}
    proc2_num_conflicts = {}
    proc2_runtime_greedy = {}
    proc2_chrom_num_greedy = {}


    while n < 1000:
        while prob <= .6:
            G = nx.gnp_random_graph(n, prob)
            #print(
            #    'Using a graph generated from G_', n, ',', prob, 'with average degree', avg_degree(G) 
            #)
            proc2_avg_deg[str(n) + ' ' + str(prob)] = avg_degree(G)

            start_time = time.process_time()
            color_p, num_conflicts = parallel_color(G, p)
            #print(
            #    'The parallelized algorithm took', time.process_time() - start_time,
            #    'seconds to run with', num_conflicts, 'conflicts and a chromatic number of',
            #    max(color_p.values())
            #    )
            proc2_runtime_parallel[str(n) + ' ' + str(prob)] = time.process_time() - start_time
            proc2_num_conflicts[str(n) + ' ' + str(prob)] = num_conflicts
            proc2_chrom_num_parallel[str(n) + ' ' + str(prob)] = max(color_p.values())

            start_time = time.process_time()
            color_g = greedy_color(G)
            #print(
            #    'Networkx greedy algorithm took', time.process_time() - start_time,
            #    'seconds to run with a chromatic number of', max(color_g.values())
            #    )
            proc2_runtime_greedy[str(n) + ' ' + str(prob)] = time.process_time() - start_time
            proc2_chrom_num_greedy = max(color_g.values())
            
            n += 100
            prob += .1

    # 4 Processors
    n = 100
    prob = .1
    p = 4
    proc4_avg_deg = {}
    proc4_runtime_parallel = {}
    proc4_chrom_num_parallel = {}
    proc4_num_conflicts = {}
    proc4_runtime_greedy = {}
    proc4_chrom_num_greedy = {}


    while n < 1000:
        while prob <= .6:
            G = nx.gnp_random_graph(n, prob)
            #print(
            #    'Using a graph generated from G_', n, ',', prob, 'with average degree', avg_degree(G) 
            #)
            proc4_avg_deg[str(n) + ' ' + str(prob)] = avg_degree(G)

            start_time = time.process_time()
            color_p, num_conflicts = parallel_color(G, p)
            #print(
            #    'The parallelized algorithm took', time.process_time() - start_time,
            #    'seconds to run with', num_conflicts, 'conflicts and a chromatic number of',
            #    max(color_p.values())
            #    )
            proc4_runtime_parallel[str(n) + ' ' + str(prob)] = time.process_time() - start_time
            proc4_num_conflicts[str(n) + ' ' + str(prob)] = num_conflicts
            proc4_chrom_num_parallel[str(n) + ' ' + str(prob)] = max(color_p.values())

            start_time = time.process_time()
            color_g = greedy_color(G)
            #print(
            #    'Networkx greedy algorithm took', time.process_time() - start_time,
            #    'seconds to run with a chromatic number of', max(color_g.values())
            #    )
            proc4_runtime_greedy[str(n) + ' ' + str(prob)] = time.process_time() - start_time
            proc4_chrom_num_greedy = max(color_g.values())
            
            n += 100
            prob += .1

    # 6 Processors
    n = 100
    prob = .1
    p = 6
    proc6_avg_deg = {}
    proc6_runtime_parallel = {}
    proc6_chrom_num_parallel = {}
    proc6_num_conflicts = {}
    proc6_runtime_greedy = {}
    proc6_chrom_num_greedy = {}


    while n < 1000:
        while prob <= .6:
            G = nx.gnp_random_graph(n, prob)
            #print(
            #    'Using a graph generated from G_', n, ',', prob, 'with average degree', avg_degree(G) 
            #)
            proc6_avg_deg[str(n) + ' ' + str(prob)] = avg_degree(G)

            start_time = time.process_time()
            color_p, num_conflicts = parallel_color(G, p)
            #print(
            #    'The parallelized algorithm took', time.process_time() - start_time,
            #    'seconds to run with', num_conflicts, 'conflicts and a chromatic number of',
            #    max(color_p.values())
            #    )
            proc6_runtime_parallel[str(n) + ' ' + str(prob)] = time.process_time() - start_time
            proc6_num_conflicts[str(n) + ' ' + str(prob)] = num_conflicts
            proc6_chrom_num_parallel[str(n) + ' ' + str(prob)] = max(color_p.values())

            start_time = time.process_time()
            color_g = greedy_color(G)
            #print(
            #    'Networkx greedy algorithm took', time.process_time() - start_time,
            #    'seconds to run with a chromatic number of', max(color_g.values())
            #    )
            proc6_runtime_greedy[str(n) + ' ' + str(prob)] = time.process_time() - start_time
            proc6_chrom_num_greedy = max(color_g.values())
            
            n += 100
            prob += .1