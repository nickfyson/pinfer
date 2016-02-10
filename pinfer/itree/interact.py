# -*- coding: utf-8 -*-
"""module docstring here"""

import networkx as nx
from networkx.exception import NetworkXNoPath

from .utils import get_inode_name
from .utils import gene_is_lost


def make_new_inode(tree, geneA, geneB):

    inode = get_inode_name(geneA, geneB)

    inode_name = get_inode_name(tree.node[geneA]['name'], tree.node[geneB]['name'])

    if inode in tree.nodes():
        raise Exception("Uh oh, this shouldn't be here!")
    tree.add_node(inode, name=inode_name, node_type='interaction', S=tree.node[geneA]['S'])

    return inode


def add_inode_parent(tree, inode):

    def walk_back_gene_pair(tree, geneA, geneB):

        parentsA = tree.predecessors(geneA)
        parentsB = tree.predecessors(geneB)

        # we check for bugs with some (theoretically) impossible cases
        if len(parentsA) == 0 and len(parentsB) == 0:
            print(nx.topological_sort(tree)[0])
            print(geneA, geneB)
            raise Exception("Neither has a parent?! This shouldn't be possible!")
        elif len(parentsA)  > 1 or len(parentsB) > 1:
            raise Exception("this shouldn't be possible!")

        # if either gene has no parents, we *must* walk back the other one
        if len(parentsA) == 0:
            return geneA, parentsB[0]
        if len(parentsB) == 0:
            return parentsA[0], geneB

        # we now know we have exactly one parent for each gene
        parentA = parentsA[0]
        parentB = parentsB[0]

        # if both genes have a parent, we must select the one most recently born
        # in the case of a tie, its abitrary which we choose
        if tree.node[parentA]['t_birth'] > tree.node[parentB]['t_birth']:
            return parentA, geneB
        else:
            return geneA, parentB

    def get_single_common_inode(tree, geneA, geneB):

        # we know that if geneA and geneB are the same, we must return the self interaction
        if geneA == geneB:
            return get_inode_name(geneA, geneB)

        # now we can check to see if the two genes have a child interaction in common
        childrenA = {n for n in tree.successors(geneA)
                     if tree.node[n]['node_type'] == 'interaction'}
        childrenB = {n for n in tree.successors(geneB)
                     if tree.node[n]['node_type'] == 'interaction'}

        common_children = childrenA.intersection(childrenB)

        if len(common_children) > 1:
            # this is a sanity check, and ought to be impossible
            print(common_children)
            raise Exception('Uh oh, you should NOT be able to have more ' +
                            'than one joint child interaction!')

        # if we don't find a common child interaction, we return false
        if len(common_children) == 0:
            return False
        # finally, if there's exactly one common child interaction, we return it
        return common_children.pop()

    genes = [n for n in tree.predecessors(inode) if tree.node[n]['node_type'] == 'gene']

    # the two parents may be the same if there is only a single predecessor
    geneA = genes[0]
    geneB = genes[-1]

    # if the genes are both the ancestral root gene, then there is no parent to find
    if geneA == geneB == nx.topological_sort(tree)[0]:
        return

    # we continue an iterative search until the parent interaction is found
    while True:

        geneA, geneB = walk_back_gene_pair(tree, geneA, geneB)

        parent_interaction = get_single_common_inode(tree, geneA, geneB)

        if parent_interaction:
            break

    tree.add_edge(parent_interaction, inode, edge_type='interaction')


def add_inode_distance(tree, iparent, inode):

    genes = [n for n in tree.predecessors(inode) if tree.node[n]['node_type'] == 'gene']

    # the two genes may be the same if this is a self interaction
    geneA = genes[0]
    geneB = genes[-1]

    parents = [n for n in tree.predecessors(iparent) if tree.node[n]['node_type'] == 'gene']

    # the two parents may be the same if this is a self interaction
    parentA = parents[0]
    parentB = parents[-1]

    if geneA == geneB:
        # there is only one path between them and parentA must also equal parentB
        assert parentA == parentB

        tree.edge[iparent][inode]['evol_dist'] = \
            nx.shortest_path_length(tree, parentA, geneA, weight='distance')
        return

    # if inode is not a self interaction, there must be two route total between
    # all (unique) parents and all genes
    path_lengths = []
    for parent in set([parentA, parentB]):
        for gene in [geneA, geneB]:
            try:
                path_lengths.append(nx.shortest_path_length(tree, parent, gene, weight='distance'))
            except NetworkXNoPath:
                pass

    if len(path_lengths) != 2:
        print(path_lengths)
        raise Exception("there aren't exactly two paths!!")

    # tree.edge[iparent][inode]['evol_dist'] = path_lengths
    tree.edge[iparent][inode]['evol_dist'] = sum(path_lengths)


def add_all_inodes(tree):
    """function to construct interaction tree, given suitably annotated gene tree"""

    from itertools import combinations_with_replacement

    # for an initial implementation, we go for the most conceptually simple approach
    # efficiency be damned!
    for (geneA, geneB) in combinations_with_replacement(tree.nodes(), 2):

        # the genes must be from the same species in order to interact
        if tree.node[geneA]['S'] != tree.node[geneB]['S']:
            continue

        # if geneA died before or at same time as geneB born, no interaction
        if tree.node[geneA]['t_death'] <= tree.node[geneB]['t_birth']:
            continue

        # if geneB died before or at same time as geneA born, no interaction
        if tree.node[geneB]['t_death'] <= tree.node[geneA]['t_birth']:
            continue

        # if either gene is lost, no interaction
        if gene_is_lost(tree, geneA) or gene_is_lost(tree, geneB):
            continue

        inode = make_new_inode(tree, geneA, geneB)

        # this new node hangs from the parent genes
        tree.add_edge(geneA, inode, distance=0.0)
        tree.add_edge(geneB, inode, distance=0.0)

    # now each inode must be connected to its parent interaction
    for inode in [n for n in tree.nodes() if tree.node[n]['node_type'] == 'interaction']:

        add_inode_parent(tree, inode)

    for s, t in [(s, t) for s, t in tree.edges()
                 if tree.edge[s][t].get('edge_type', '') == 'interaction']:

        add_inode_distance(tree, s, t)

    # we don't want the gTree nodes actually remaining as part of the iTree
    tree.remove_nodes_from([n for n in tree.nodes() if tree.node[n]['node_type'] == 'gene'])

    # instead of inscrutable numbers as the nodes,
    # we relabel using the 'name' property
    nx.relabel_nodes(tree,
                     {n: tree.node[n]['name'] for n in tree.nodes()},
                     copy=False)

    return tree
