# -*- coding: utf-8 -*-
"""module docstring here"""


import networkx as nx


def analyse_tree(tree, update_repeats=20, verbose=True):
    """perform VB inference on networkx representation of interaction tree"""

    from bayespy.nodes import Categorical, Mixture, Gate

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

    inf_engine = VB(*variables.values())

    # any node with the property 'observed' has its value fixed
    for node in tree.nodes():
        if 'observed' in tree.node[node]:
            variables[node].observe(tree.node[node]['observed'])

    inf_engine.update(repeat=update_repeats, verbose=verbose)

    for node, variable in variables.items():
        tree.node[node]['posterior'] = variable.pdf(1)


def analyse_pymc(tree, samples=1000, burns=500):
    """docstring for analyse_pymc"""
    
    import numpy as np
    import pymc as pm

    variables = {}

    for node in nx.topological_sort(tree):

        # if there are no incoming edges, this node is the initial state
        if tree.in_degree(node) == 0:
            prior = tree.node[node]['prior']
            
            variables[node] = pm.Categorical(node, [1 - prior, prior])

        for s, t in tree.in_edges(node):

            p_transitions = tree.edge[s][t]['p_transitions']
            
            variables['t_%s' % t] = pm.Index('t_%s' % t, p_transitions, variables[s])
            
            if 'observed' in tree.node[t]:
                observed = True
                value    = int(tree.node[t]['observed'])
            else:
                observed = None
                value    = None
            
            variables[t] = pm.Categorical(t, variables['t_%s' % t],
                                          observed=observed, value=value)
    
    model = pm.MCMC(list(variables.values()))

    model.sample(10000, 5000)

    for node in tree.nodes():
        tree.node[node]['posterior'] = np.mean(model.trace(node)[:])

    model = None

    return tree
