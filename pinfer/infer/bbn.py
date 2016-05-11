# -*- coding: utf-8 -*-
"""module docstring here"""

import numpy as np
import pickle
from hashlib import sha224


def get_function_name(string):

    return 'func_' + sha224(pickle.dumps(string)).hexdigest()


def get_variable_name(argument):

    if isinstance(argument, list):
        return ['var_' + sha224(pickle.dumps(string)).hexdigest() for string in argument]
    else:
        return 'var_' + sha224(pickle.dumps(argument)).hexdigest()


def create_node_functions(tree):

    f_lines = []

    for node in tree.nodes():

        fname = get_function_name(node)

        arguments = sorted(tree.pred[node]) + [node]

        try:
            p_array = tree.node[node]['CPT']
        except KeyError:
            p_array = tree.node[node]['prior']

        f_lines.append("def %s(%s):" % (fname, ','.join(get_variable_name(arguments))))
        # f_lines.append("    '''Cancer'''")
        f_lines.append("    table = dict()")
        for indices, value in np.ndenumerate(p_array):
            f_lines.append("    table['%s'] = %f" % (''.join([str(x) for x in indices]), value))

        f_lines.append("    key = ''")
        for argument in arguments:
            f_lines.append("    key = key + '1' if %s else key + '0'" %
                           (get_variable_name(argument)))
        f_lines.append("    return table[key]")
        f_lines.append("")

    functions_string = '\n'.join(f_lines)

    functions = {}

    exec(functions_string, functions)

    return [functions[key] for key in functions.keys() if '__' not in key]


def analyse_bbn(tree):
    """

    """

    # from .. import bayesian
    # from bayesian.bbn import build_bbn
    from bayesian.factor_graph import build_graph

    functions = create_node_functions(tree)

    # g = build_bbn(functions)
    g = build_graph(functions)

    observations = {}
    for node in [n for n in tree.nodes() if 'observation' in tree.node[n]]:
        if (tree.node[node]['observation'] == np.array([0.0, 1.0])).all():
            observations[get_variable_name(node)] = 1
        elif (tree.node[node]['observation'] == np.array([1.0, 0.0])).all():
            observations[get_variable_name(node)] = 0
        else:
            raise Exception('Have to have binary observations!')

    results = g.query(**observations)

    for node in tree.nodes():
        tree.node[node]['belief'] = np.zeros(2)

        tree.node[node]['belief'][0] = results[(get_variable_name(node), False)]
        tree.node[node]['belief'][1] = results[(get_variable_name(node), True)]

    return tree
