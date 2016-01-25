# -*- coding: utf-8 -*-
"""module docstring here"""

from copy import deepcopy
import networkx as nx

from .utils import get_inode_name, gene_is_lost
from .evol_time import label_birth_death


def get_fellow_extants(iTree, gene):
    """returns a list of all interaction partners of the specified gene"""

    if gene_is_lost(iTree, gene):
        return []

    time = iTree.node[gene]['t_birth']

    extants = []

    for n in iTree.nodes():
        if iTree.node[n]['node_type'] != 'gene':
            pass
        elif gene_is_lost(iTree, n):
            pass
        elif iTree.node[n]['t_birth'] > time:
            pass
        elif iTree.node[n]['t_death'] <= time:
            pass
        elif iTree.node[n]['S'] != iTree.node[gene]['S']:
            pass
        else:
            extants.append(n)

    return extants


def get_parent_interaction(iTree, interaction):
    """return the appropriate parent for the given interaction"""

    # the parent gene(s) are attached as predecessors in the graph

    predecessors = iTree.predecessors(interaction)

    geneA = predecessors[0]
    try:
        geneB = predecessors[1]
    except IndexError:  # if there isn't a 2nd predecessor, we know it's a self interaction
        geneB = predecessors[0]

    try:
        if iTree.node[geneA]['t_birth'] > iTree.node[geneB]['t_birth']:
            ancestorA = iTree.predecessors(geneA)[0]
            ancestorB = geneB
        else:
            ancestorA = geneA
            ancestorB = iTree.predecessors(geneB)[0]
    except IndexError:
        # if there are no ancestors, then we're looking at the ancenstral interaction
        return None, None

    while True:

        # if we have found the same ancestor for both, we know that
        # the parent interaction must be the self-interaction
        if ancestorA == ancestorB:
            parent_interaction = get_inode_name(ancestorA, ancestorB)
            break

        # find the common interaction child of these two parent nodes
        childrenA = set([n for n in iTree.successors(ancestorA)
                         if iTree.node[n]['node_type'] == 'interaction'])
        childrenB = set([n for n in iTree.successors(ancestorB)
                         if iTree.node[n]['node_type'] == 'interaction'])

        # as a sanity check, it should be impossible to have *more* than 1 common child
        if len(childrenA.intersection(childrenB)) > 1:
            raise Exception('Uh oh! Too many children in common!')
        # if there is a common node, we have a winner!
        if childrenA.intersection(childrenB):
            parent_interaction = childrenA.intersection(childrenB).pop()
            break

        # if not, we replace the most recently deceased ancestor node with *its* ancestor gene
        if iTree.node[ancestorA]['t_birth'] > iTree.node[ancestorB]['t_birth']:
            ancestorA = iTree.predecessors(ancestorA)[0]
            ancestorB = ancestorB
        else:
            ancestorA = ancestorA
            ancestorB = iTree.predecessors(ancestorB)[0]

    # we can also return the edge length here...
    evolution_of_A = nx.shortest_path_length(iTree, source=ancestorA,
                                             target=geneA, weight='distance')
    evolution_of_B = nx.shortest_path_length(iTree, source=ancestorB,
                                             target=geneB, weight='distance')
    evol_dist    = evolution_of_A + evolution_of_B

    return parent_interaction, evol_dist


def build_itree(gTree):
    """function to construct interaction tree, given suitably annotated gene tree"""

    # we initialise as a copy of the gTree, so the new tree can be built on the existing structure
    iTree = deepcopy(gTree)
    iTree.graph['name'] = 'iTree'

    # all existing nodes are annotated as genes
    for node in iTree.nodes():
        iTree.node[node]['node_type'] = 'gene'

    label_birth_death(iTree)

    # first we build a list of nodes ordered by their birth time
    nodes_t_births = [(n, iTree.node[n]['t_birth']) for n in iTree.nodes()]
    ordered_nodes  = [n for n, d in sorted(nodes_t_births, key=lambda x:x[1])]

    for gene in ordered_nodes:

        for fellow in get_fellow_extants(iTree, gene):

            # if doesn't already exist, add an interaction between node and gene
            if get_inode_name(gene, fellow) not in iTree.nodes():

                new_interaction = get_inode_name(gene, fellow)
                new_int_name    = get_inode_name(gTree.node[gene]['name'],
                                                 gTree.node[fellow]['name'])
                iTree.add_node(new_interaction,
                               node_type='interaction',
                               S=iTree.node[gene]['S'],
                               name=new_int_name
                               )

                iTree.add_edge(gene, new_interaction)
                iTree.add_edge(fellow, new_interaction)

                parent_interaction, evol_dist = get_parent_interaction(iTree, new_interaction)

                if parent_interaction:
                    iTree.add_edge(parent_interaction, new_interaction, evol_dist=evol_dist)

    # we can use the remaining gTree nodes to mark all the extant interactions
    # we first list all the extant gene nodes in the present time - ie. no children
    extant_gnodes = set([n for n in gTree.nodes() if not gTree.successors(n)])

    for inode in [n for n in iTree.nodes() if iTree.node[n]['node_type'] == 'interaction']:

        parent_genes = [n for n in iTree.predecessors(inode)
                        if iTree.node[n]['node_type'] == 'gene']
        # if the parent gene(s) is / are both extant, then so is the interaction
        if set(parent_genes).issubset(extant_gnodes):
            iTree.node[inode]['extant'] = True

    # we don't want the gTree nodes actually remaining as part of the iTree
    iTree.remove_nodes_from([n for n in iTree.nodes() if iTree.node[n]['node_type'] == 'gene'])

    # instead of inscrutable numbers as the nodes,
    # we relabel using the 'name' property
    nx.relabel_nodes(iTree,
                     {n: iTree.node[n]['name'] for n in iTree.nodes()},
                     copy=False)

    return iTree
