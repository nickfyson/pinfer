# -*- coding: utf-8 -*-
"""module docstring here"""

import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx


def vis_tree(tree, fig=None, pos=None, node_cat_string=None, node_color_string=None,
             node_size=None, layout='dot'):
    """visualise the tree with colouring according to property 'node_cat_string'"""

    if not fig:
        plt.figure(figsize=(24, 8))

    if not pos:
        pos = nx.graphviz_layout(tree, prog=layout)

    node_colors = []

    if node_cat_string:

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

    if node_color_string:

        node_colors = []

        for n in tree.nodes():
            node_colors.append(tree.node[n].get(node_color_string, (0.8, 0.8, 0.8)))

    if not node_size:
        node_size = float(6e4) / len(tree.nodes())

    if node_size > 2000:
        with_labels = True
        font_size   = node_size / 400.0
    else:
        with_labels = False
        font_size   = node_size / 400.0

    if not node_colors:
        node_colors = [(0.8, 0.8, 0.8) for n in tree.nodes()]

    nx.draw(tree, pos, arrows=False,
            with_labels=with_labels, font_size=font_size,
            edge_color=[tree.edge[s][t].get('color', '#000000') for s, t in tree.edges()],
            width=[tree.edge[s][t].get('width', 1.0) for s, t in tree.edges()],
            node_color=node_colors,
            linewidths=0.0,
            alpha=0.8,
            node_size=node_size)
