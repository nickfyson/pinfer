# -*- coding: utf-8 -*-
"""module docstring here"""

import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx


def vis_tree(tree, fig=None, pos=None, node_cat_string=None,
             pallete_name='rainbow', node_size=None, layout='dot'):
    """visualise the tree with colouring according to property 'node_cat_string'"""

    if not fig:
        plt.figure(figsize=(20, 6))

    if not pos:
        pos = nx.graphviz_layout(tree, prog=layout)

    if node_cat_string:

        all_categories = set([tree.node[n][node_cat_string] for n in tree.nodes()
                              if tree.node[n].get(node_cat_string, None)])

        all_categories = sorted(list(all_categories))

        pallette = sns.color_palette(pallete_name, len(all_categories))

        color_dict = {cat: pallette[i] for i, cat in enumerate(all_categories)}

        node_colors = []
        for node in tree.nodes():

            category = tree.node[node].get(node_cat_string, '')

            node_colors.append(color_dict.get(category, (0.8, 0.8, 0.8)))

    if not node_size:
        node_size = float(6e4) / len(tree.nodes())

    if node_size > 2000:
        with_labels = True
        font_size   = node_size / 400.0
    else:
        with_labels = False
        font_size   = node_size / 400.0

    nx.draw(tree, pos, arrows=False,
            with_labels=with_labels, font_size=font_size,
            node_color=node_colors,
            alpha=0.8,
            linewidths=0.0, node_size=node_size)
