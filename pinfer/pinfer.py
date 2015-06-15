# -*- coding: utf-8 -*-

import networkx as nx


def analyse_tree(tree, update_repeats=20, verbose=True):
    """perform VB inference on networkx representation of interaction tree"""

    from bayespy.nodes import Categorical, Mixture

    variables = {}

    for node in nx.topological_sort(tree):

        # if there are no incoming edges, this node is the initial state
        if tree.in_degree(node) == 0:
            prior = tree.node[node]['prior']

            variables[node] = Categorical([1 - prior, prior])

        for origin, target in tree.in_edges(node):

            p_transitions = tree.edge[origin][target]['p_transitions']

            variables[target] = Mixture(variables[origin],
                                        Categorical, p_transitions)

    from bayespy.inference import VB

    Q = VB(*variables.values())

    # any node with the property 'observed' has its value fixed
    for node in tree.nodes():
        if 'observed' in tree.node[node]:
            variables[node].observe(tree.node[node]['observed'])

    Q.update(repeat=update_repeats, verbose=verbose)

    for node, variable in variables.items():
        tree.node[node]['posterior'] = variable.pdf(1)
