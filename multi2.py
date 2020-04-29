from multiprocessing import Pool, Manager, Process
import networkx as nx
import numpy as np
import itertools

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


def nodes_iter(Vp, p):
    ''' Given p chunks of nodes, order the vertices into an array, blocks, such that blocks[i]
    indicates vertices that will be colored in parallel at iteration i '''
    min_blen = min([len(block) for block in Vp])
    blocks = []

    for i in range(min_blen):
        index_dic = {}
        for j in range(p):
            if i not in index_dic:
                index_dic[i] = [Vp[j][i]]
            else:
                index_dic[i] += [Vp[j][i]]
        blocks.append(index_dic)

    index_dic = {}
    for chunk in Vp:
        if len(chunk) > min_blen:
            i = min_blen
            if i not in index_dic:
                index_dic[i] = [chunk[i]]
            else:
                index_dic[i] += [chunk[i]]
    blocks.append(index_dic)
    return blocks


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
    blocks = nodes_iter(Vp, min_blen)

    coloring = {}

    for i in range(len(blocks)):
        p1 = Pool(p)
        for u in blocks[i][i]:
            Gnew = induce_subgraph(G, u, blocks[i][i])
            (u, color) = (p1.apply_async(func = color_node, args = (Gnew, u, coloring,))).get()
            coloring[u] = color
    
    # Phase 2
    A = []

    for i in range(len(blocks)):
        for v in blocks[i][i]:
            for u in G.neighbors(v):
                if u in blocks[i][i]:
                    if coloring[u] == coloring[v]:
                        if min(u, v) not in A:
                            A.append(min(u, v))

    for v in A:
        color_node(G, v, coloring)
    return coloring



if __name__ == "__main__":
    G = nx.gnp_random_graph(5, .5)
    p = 2
    print(parallel_color(G, p), G.edges)