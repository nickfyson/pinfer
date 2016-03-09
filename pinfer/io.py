# -*- coding: utf-8 -*-
"""module docstring here"""

import re

import networkx as nx

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

    graph = nx.DiGraph()

    for node in tree.nodes():
        new_node = node_translator[node]

        properties = {'name': str(node)}
        for match in re.findall(r'[^:]*\=[^:]*', node.comment):
            properties[match.split('=')[0]] = match.split('=')[1]

        graph.add_node(new_node, **properties)

    for source, target in tree.edges():
        new_source = node_translator[source]
        new_target = node_translator[target]
        graph.add_edge(new_source, new_target,
                       distance=source.distance(target),
                       **tree.edge[source][target])

    for s, t in graph.edges():
        graph.edge[s][t].pop('weight')

    # follow convention by renaming the root node to 'X0'
    root                     = nx.topological_sort(graph)[0]
    graph.node[root]['name'] = 'X0'

    # rename lost genes, so all nodes have unique names
    for n in [n for n in graph.nodes() if 'lost' in graph.node[n]['name'].lower()]:
        graph.node[n]['name'] = n + graph.node[n]['name']

    # build dictionary to replace node objects with the name of each node
    new_node_names = {n: graph.node[n]['name'] for n in graph.nodes()}

    nx.relabel_nodes(graph, new_node_names, copy=False)

    return graph
