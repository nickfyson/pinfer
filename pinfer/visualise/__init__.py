# -*- coding: utf-8 -*-
"""module docstring here"""

import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx


def vis_tree(tree, fig=None, node_cat_string=None,
             node_size=None, layout='dot', alpha=0.8):
    """visualise the tree with colouring according to property 'node_cat_string'"""

    if not fig:
        plt.figure(figsize=(24, 8))

    pos = tree.graph.get('pos', None)
    if not pos:
        pos = nx.graphviz_layout(tree, prog=layout)

    if node_cat_string:

        node_colors = []

        all_categories = set([tree.node[n][node_cat_string] for n in tree.nodes()
                              if tree.node[n].get(node_cat_string, None)])

        all_categories = sorted(list(all_categories))

        pallette = sns.cubehelix_palette(n_colors=len(all_categories), start=1.1,
                                         dark=0.4, light=0.8, rot=2.5, hue=0.9, gamma=1.0)

        color_dict = {cat: pallette[i] for i, cat in enumerate(all_categories)}

        node_colors = []
        for node in tree.nodes():

            category = tree.node[node].get(node_cat_string, '')

            node_colors.append(color_dict.get(category, (0.8, 0.8, 0.8)))
    else:

        node_colors = [tree.node[n].get('color', (0.8, 0.8, 0.8)) for n in tree.nodes()]

    linewidths  = [tree.node[n].get('linewidth', 0.0) for n in tree.nodes()]

    default_node_size = float(6e4) / len(tree.nodes())

    node_sizes = [tree.node[n].get('size', default_node_size) for n in tree.nodes()]

    edge_colors = [tree.edge[s][t].get('color', '#000000') for s, t in tree.edges()]
    edge_widths = [tree.edge[s][t].get('width', 0.5) for s, t in tree.edges()]

    if default_node_size > 2000:
        with_labels = True
        font_size   = default_node_size / 400.0
    else:
        with_labels = False
        font_size   = default_node_size / 400.0

    nx.draw(tree, pos, arrows=False,
            with_labels=with_labels, font_size=font_size,
            edge_color=edge_colors,
            width=edge_widths,
            node_color=node_colors,
            linewidths=linewidths,
            alpha=alpha,
            node_size=node_sizes)
