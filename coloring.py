import networkx as nx
import numpy as np
import multiprocessing
from multiprocessing import Pool
import itertools

# Initialize graph from an np matrix (np.matrix())
def CreateGraph(A):
	G = nx.Graph()
	G= nx.from_numpy_matrix(A)
	return G
    
def strategy_largest_first(G, colors):
    """Returns a list of the nodes of ``G`` in decreasing order by
    degree.

    ``G`` is a NetworkX graph. ``colors`` is ignored.

    """
    return sorted(G, key=G.degree, reverse=True)

STRATEGIES = {'largest_first': strategy_largest_first}

def greedy_color(G, strategy='largest_first', interchange=False):
    if len(G) == 0:
        return {}
    # Determine the strategy provided by the caller.
    strategy = STRATEGIES.get(strategy, strategy)
    if not callable(strategy):
        raise nx.NetworkXError('strategy must be callable or a valid string. '
                               '{0} not valid.'.format(strategy))
    # Perform some validation on the arguments before executing any
    # strategy functions.
    if interchange:
        if strategy is strategy_independent_set:
            msg = 'interchange cannot be used with independent_set'
            raise nx.NetworkXPointlessConcept(msg)
        if strategy is strategy_saturation_largest_first:
            msg = ('interchange cannot be used with'
                   ' saturation_largest_first')
            raise nx.NetworkXPointlessConcept(msg)
    colors = {}
    nodes = strategy(G, colors)
    if interchange:
        return _interchange.greedy_coloring_with_interchange(G, nodes)
    for u in nodes:
        # Set to keep track of colors of neighbours
        neighbour_colors = {colors[v] for v in G[u] if v in colors}
        # Find the first unused color.
        for color in itertools.count():
            if color not in neighbour_colors:
                break
        # Assign the new color to the current node.
        colors[u] = color
    return colors


def Coloring1(G, p):
    # Partition vertices into p blocks
    Vp = [list(G.nodes)[i:i + p] for i in range(0, len(G.nodes), p)]
    blocks = [G.subgraph(Vp[i]) for i in range(p)]
    coloring = {}
    for i in range(p):
        coloring[blocks[i]] = greedy_color(blocks[i])
    
            
    
if __name__ == "__main__":
    A = np.matrix([[1,1, 1],[1,1, 1], [1, 1, 1]])
    G = CreateGraph(A)
    print(Coloring1(G, 2))
