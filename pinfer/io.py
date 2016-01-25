# -*- coding: utf-8 -*-
"""module docstring here"""

import re

from networkx import DiGraph

from Bio.Phylo import read, to_networkx


def load_notung_nhx(filename):
    """load reconciled gene tree from NHX formatted file

    returns networkx graph object
    strips information from the comment field and converts into node properties"""

    with open(filename, 'r') as f:
        tree = read(f, format='newick')

    tree.rooted = True

    tree = to_networkx(tree)

    node_translator = {}
    for node in tree.nodes():
        node_translator[node] = str(len(node_translator))

    graph = DiGraph()

    for node in tree.nodes():
        new_node = node_translator[node]

        properties = {'name': str(node)}
        for match in re.findall(r'[^:]*\=[^:]*', node.comment):
            properties[match.split('=')[0]] = match.split('=')[1]

        graph.add_node(new_node, **properties)

    for source, target in tree.edges():
        new_source = node_translator[source]
        new_target = node_translator[target]
        graph.add_edge(new_source, new_target, **tree.edge[source][target])

    for s, t in graph.edges():
        graph.edge[s][t]['distance'] = graph.edge[s][t].pop('weight')

    return graph
