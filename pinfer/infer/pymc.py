
import networkx as nx


def analyse_pymc(tree, samples=1000, burns=500):
    """docstring for analyse_pymc"""

    import numpy as np
    import pymc as pm

    nodes = {}
    trans = {}

    for node in nx.topological_sort(tree):

        # if there are no incoming edges, this node is the initial state
        if tree.in_degree(node) == 0:
            prior = tree.node[node]['prior']

            nodes[node] = pm.Categorical(node, [1 - prior, prior])

        for s, t in tree.in_edges(node):

            p_transitions = tree.edge[s][t]['p_transitions']

            trans[t] = pm.Index('tran_%s' % t, p_transitions, nodes[s])

            if 'observed' in tree.node[t]:
                nodes[t] = pm.Categorical(t, trans[t],
                                          observed=True, value=tree.node[t]['observed'])
            else:
                nodes[t] = pm.Categorical(t, trans[t])

    model = pm.MCMC(list(nodes.values()) + list(trans.values()))

    model.sample(10000, 5000)

    for node in tree.nodes():
        tree.node[node]['posterior'] = np.mean(model.trace(node)[:])

    model = None

    return tree
