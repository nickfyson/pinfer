# -*- coding: utf-8 -*-
"""module docstring here"""

from copy import deepcopy
import networkx as nx


def get_inode_name(geneA, geneB):
    # the name of interaction nodes is a concatenation of the two gene names
    # crucially, these are always sorted so the order in which genes are passed is irrelevant
    return '%s_%s' % tuple(sorted((geneA, geneB)))


def get_ordered_nodes(tree, include_leaves=True):
    """return list of all nodes, ordered by distance from tree
    optionally, also return distances as a dictionary"""

    root = nx.topological_sort(tree)[0]

    nodes_dists = []
    for node in tree.nodes():
        if tree.out_degree(node) == 0 and not include_leaves:
            pass
        else:
            nodes_dists.append((nx.shortest_path_length(tree, root, node, weight='length'), node))

    ordered_nodes = [node for (dist, node) in sorted(nodes_dists)]

    distances = {node: dist for (dist, node) in sorted(nodes_dists)}

    return ordered_nodes, distances


def add_normalised_edge_lengths(tree):
    """annotate all edges such that total depth of tree is normalised to one within each species"""

    def relabel_weights(tree, effective_root, remaining_length=1.0):
        
        for child in tree.edge[effective_root]:
            
            max_dist = tree.edge[effective_root][child]['weight'] + \
                max([
                    d for n, d in
                    nx.shortest_path_length(tree, source=child, weight='weight').items()
                ])
            
            original_length = tree.edge[effective_root][child]['weight']
            
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


def label_death_times(tree):
    """add birth and death time properties to each node"""

    add_normalised_edge_lengths(tree)
    
    root = [n for n in tree.nodes() if not tree.predecessors(n)][0]

    for n, d in nx.shortest_path_length(tree, root, weight='length').items():
        tree.node[n]['t_death'] = d + 1.0

    for n in tree.nodes():
        try:
            tree.node[n]['t_birth'] = tree.node[tree.predecessors(n)[0]]['t_death']
        except IndexError:
            tree.node[n]['t_birth'] = 0.0
    
    return
    

def get_extants(iTree, gene):
    """returns a list of all interaction partners of the specified gene"""

    if _is_lost(iTree, gene):
        return []

    time = iTree.node[gene]['t_birth']
    
    extants = []
    
    for n in iTree.nodes():
        if iTree.node[n]['node_type'] != 'gene':
            pass
        elif _is_lost(iTree, n):
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
    
    # geneA, geneB = interaction.split('_')
    predecessors = iTree.predecessors(interaction)
    geneA = predecessors[0]
    try:
        geneB = predecessors[1]
    except IndexError:
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
    evolution_of_A = nx.shortest_path_length(iTree, source=ancestorA, target=geneA, weight='weight')
    evolution_of_B = nx.shortest_path_length(iTree, source=ancestorB, target=geneB, weight='weight')
    evol_dist    = evolution_of_A + evolution_of_B

    return parent_interaction, evol_dist
    

def _is_lost(iTree, gene):
    """simple function to determine whether the gene has been lost"""
    return 'lost' in iTree.node[gene]['name'].lower()


def build_itree(gTree):
    """function to construction interaction tree, given suitably annotated gene tree"""
    # we initialise as a copy of the gTree, so the new tree can be built on the existing structure
    
    iTree = deepcopy(gTree)
    iTree.graph['name'] = 'iTree'

    # all existing nodes are annotated as genes
    for node in iTree.nodes():
        iTree.node[node]['node_type'] = 'gene'
 
    label_death_times(iTree)

    # first we build a list of nodes ordered by the distance from the root
    # ordered_nodes = get_ordered_nodes(iTree)
    
    ordered_nodes = [(n, iTree.node[n]['t_birth']) for n in iTree.nodes()]

    ordered_nodes = [n for n, d in sorted(ordered_nodes, key=lambda x:x[1])]

    for gene in ordered_nodes:

        for extant in get_extants(iTree, gene):
            
            # if doesn't already exist, add an interaction between node and gene
            if get_inode_name(gene, extant) not in iTree.nodes():
                
                new_interaction = get_inode_name(gene, extant)
                new_int_name    = get_inode_name(gTree.node[gene]['name'],
                                                 gTree.node[extant]['name'])
                iTree.add_node(new_interaction,
                               node_type='interaction',
                               S=iTree.node[gene]['S'],
                               name=new_int_name
                               )

                iTree.add_edge(gene, new_interaction)
                iTree.add_edge(extant, new_interaction)

                parent_interaction, evol_dist = get_parent_interaction(iTree, new_interaction)

                if parent_interaction:
                    iTree.add_edge(parent_interaction, new_interaction, evol_dist=evol_dist)
    
    # we can use the remaining gTree nodes to mark all the extant interactions
    extant_gnodes = set([n for n in gTree.nodes() if not gTree.successors(n)])

    for inode in [n for n in iTree.nodes() if iTree.node[n]['node_type'] == 'interaction']:

        parent_genes = [n for n in iTree.predecessors(inode)
                        if iTree.node[n]['node_type'] == 'gene']
               
        if set(parent_genes).issubset(extant_gnodes):
            iTree.node[inode]['extant'] = True
    
    # we don't want the gTree nodes actually remaining as part of the iTree
    iTree.remove_nodes_from([n for n in iTree.nodes() if iTree.node[n]['node_type'] == 'gene'])
    
    return iTree
