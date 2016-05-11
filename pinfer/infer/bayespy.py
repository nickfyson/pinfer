
import networkx as nx
import numpy as np
from bayespy.nodes import Categorical, MultiMixture


def analyse_bayespy(tree, update_repeats=100, verbose=False):
    """perform VB inference on networkx representation of interaction tree"""

    if nx.cycle_basis(tree.to_undirected()):
        raise Exception('BayesPy can (currently) only be used on polytrees! ' +
                        '(no cycles when undirected)')

    variables = {}

    for node in nx.topological_sort(tree):

        # if there are no incoming edges, this node is the initial state
        if tree.in_degree(node) == 0:

            variables[node] = Categorical(tree.node[node]['prior'])

        else:

            parents = tuple(variables[n] for n in sorted(tree.predecessors(node)))

            variables[node] = MultiMixture(parents, Categorical, tree.node[node]['CPT'].tolist())

    for node in [n for n in tree.nodes() if 'observation' in tree.node[n]]:
        if (tree.node[node]['observation'] == np.array([0.0, 1.0])).all():
            variables[node].observe(1)
        elif (tree.node[node]['observation'] == np.array([1.0, 0.0])).all():
            variables[node].observe(0)
        else:
            raise Exception('Have to have binary observations!')

    from bayespy.inference import VB
    inf_engine = VB(*list(variables.values()))
    inf_engine.update(repeat=update_repeats, verbose=verbose)

    for node, variable in variables.items():
        tree.node[node]['belief'] = [variable.pdf(0), variable.pdf(1)]

    return tree


    variables = {}

    for node in nx.topological_sort(tree):

        # if there are no incoming edges, this node is the initial state
        if tree.in_degree(node) == 0:
            prior = tree.node[node]['prior'].tolist()

            variables[node] = Categorical(prior)

        for s, t in tree.in_edges(node):

            variables[s]

            p_transitions = tree.node[t]['CPT'].tolist()

            variables[target] = Mixture(variables[origin],
                                        Categorical, p_transitions)

    from bayespy.inference import VB

    inf_engine = VB(*variables.values())

    # any node with the property 'observed' has its value fixed
    for node in tree.nodes():
        if 'observed' in tree.node[node]:
            variables[node].observe(tree.node[node]['observed'])

    inf_engine.update(repeat=update_repeats, verbose=verbose)

    for node, variable in variables.items():
        tree.node[node]['belief'] = [variable.pdf(0), variable.pdf(1)]
