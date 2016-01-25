
import networkx as nx


def add_normalised_edge_lengths(tree):
    """annotate all edges such that total depth of tree is normalised to one within each species"""

    def relabel_weights(tree, effective_root, remaining_length=1.0):

        for child in tree.edge[effective_root]:

            max_dist = tree.edge[effective_root][child]['distance'] + \
                max([
                    d for n, d in
                    nx.shortest_path_length(tree, source=child, weight='distance').items()
                ])

            original_length = tree.edge[effective_root][child]['distance']

            new_length      = remaining_length * (original_length / max_dist)

            # label outgoing edges with the appropriate fraction of the remaining length
            tree.edge[effective_root][child]['length'] = new_length

            # call relabel_weights with each child as effective root, and correct remaining length
            relabel_weights(tree, child, remaining_length=(remaining_length - new_length))

        return

    all_species = set([tree.node[n]['S'] for n in tree.nodes()])

    for i, species in enumerate(all_species):

        nodes = [n for n in tree.nodes() if tree.node[n]['S'] == species]

        # add parents to the set of nodes
        parents = set()
        for node in nodes:
            if tree.predecessors(node):
                parents.add(tree.predecessors(node)[0])

        subTree = nx.subgraph(tree, nodes + list(parents))

        for root in [n for n in subTree.nodes() if not subTree.predecessors(n)]:
            relabel_weights(subTree, root)

        for s, t in subTree.edges():
            tree.edge[s][t]['length'] = subTree.edge[s][t]['length']

    return


def label_birth_death(tree):
    """add birth and death time properties to each node"""

    root = [n for n in tree.nodes() if not tree.predecessors(n)][0]

    for n, d in nx.shortest_path_length(tree, root, weight='length').items():
        tree.node[n]['t_death'] = d + 1.0

    for n in tree.nodes():
        try:
            tree.node[n]['t_birth'] = tree.node[tree.predecessors(n)[0]]['t_death']
        except IndexError:
            tree.node[n]['t_birth'] = 0.0

    return
