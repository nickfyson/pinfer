# -*- coding: utf-8 -*-
"""module docstring here"""

# import sys
# import networkx as nx
import numpy as np
# from copy import deepcopy

def create_node_functions(tree):

    f_lines = []

    for node in tree.nodes():

        fname = node
        arguments = sorted(tree.pred[node]) + [node]
        try:
            p_array = tree.node[node]['CPT']
        except KeyError:
            p_array = tree.node[node]['prior']

        f_lines.append("def f%s(%s):" % (fname, ','.join(arguments)))
        # f_lines.append("    '''Cancer'''")
        f_lines.append("    table = dict()")
        for indices, value in np.ndenumerate(p_array):
            f_lines.append("    table['%s'] = %f" % (''.join([str(x) for x in indices]), value))

        f_lines.append("    key = ''")
        for argument in arguments:
            f_lines.append("    key = key + '1' if %s else key + '0'" % (argument))
        f_lines.append("    return table[key]")
        f_lines.append("")

    functions_string = '\n'.join(f_lines)

    functions = {}

    exec(functions_string, functions)

    return [functions[key] for key in functions.keys() if '__' not in key]

def analyse_bbn(tree):
    """

    """

    from bayesian.bbn import build_bbn

    functions = create_node_functions(tree)

    g = build_bbn(functions)

    observations = {}
    for node in [n for n in tree.nodes() if 'observation' in tree.node[n]]:
        if (tree.node[node]['observation'] == np.array([0.0, 1.0])).all():
            observations[node] = 1
        elif (tree.node[node]['observation'] == np.array([1.0, 0.0])).all():
            observations[node] = 0
        else:
            raise Exception('Have to have binary observations!')

    results = g.query(**observations)

    for node in tree.nodes():
        tree.node[node]['belief'] = np.zeros(2)

    for (node, state), belief in results.items():
        tree.node[node]['belief'][state] = belief

    return tree
