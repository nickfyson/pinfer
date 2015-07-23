# -*- coding: utf-8 -*-
"""module docstring here"""

from copy import deepcopy
import networkx as nx
from collections import defaultdict


def _inode_name(geneA, geneB):
    # the name of interaction nodes is a concatenation of the two gene names
    # crucially, these are always sorted so the order in which genes are passed is irrelevant
    return '%s_%s' % tuple(sorted((geneA, geneB)))


def get_ordered_nodes(tree, include_leaves=False):
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


def normalise_edge_lengths(tree):
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


def find_extant_interacting_genes(iTree, gene):
    """return a list of all genes that are both extant and have interactions with passed gene"""
    print('gene = ', gene)
    # first we retrieve a list of all the interaction nodes connected to the gene
    interactions = [n for n in iTree.successors(gene)
                    if iTree.node[n]['node_type'] == 'interaction']
    print('interactions = ', interactions)
    # now we find all the parent nodes of these interactions
    potentials = []
    for i in interactions:
        potentials += [n for n in iTree.predecessors(i) if iTree.node[n]['node_type'] == 'gene']
    print('potentials = ', potentials)
    # now we filter the potentials ensuring they are extant
    potentials = [n for n in potentials if iTree.node[n]['extant'] and iTree.node[n]['D'] == 'Y']
    print('potentials = ', potentials)
    # finally, we make these unique and return the set
    print('unique potentials = ', set(potentials))
    return set(potentials)


def _is_lost(iTree, gene):
    """simple function to determine whether the gene has been lost"""
    return 'lost' in iTree.node[gene]['name'].lower()


def build_itree(gTree):
    """function to construction interaction tree, given suitably annotated gene tree"""
    # we initialise as a copy of the gTree, so the new tree can be built on the existing structure
    
    iTree = deepcopy(gTree)
    iTree.graph['name'] = 'iTree'
    iTree.graph['extants'] = defaultdict(set)

    # all existing nodes are annotated as genes
    for node in iTree.nodes():
        iTree.node[node]['node_type'] = 'gene'

    normalise_edge_lengths(iTree)

    # first we build a list of nodes ordered by the distance from the root
    ordered_nodes, distances = get_ordered_nodes(iTree)

    root = ordered_nodes[0]
    # we initialise with the ancestral self-interaction
    interaction = _inode_name(root, root)
    iTree.add_node(interaction, node_type='interaction', S=iTree.node[root]['S'])
    # we begin with a single extant gene, and add it to the relevant species set
    iTree.graph['extants'][iTree.node[root]['S']].add(root)

    # we now sweep through the original gene tree,
    # moving the effective time forward till we hit each branching point
    # in each iteration we remove the node, and add interactions from its children
    for dying_gene in ordered_nodes:

        # if _is_lost(iTree, dying_gene):
            # continue
        
        # newly born genes are the child genes of the dying gene
        born_genes = [n for n in iTree.edge[dying_gene] if iTree.node[n]['node_type'] == 'gene']
        
        # duplication and non-duplication events are handled differently
        if iTree.node[dying_gene]['D'] == 'Y':

            # the parent node has branched, and is no longer extant
            iTree.graph['extants'][iTree.node[dying_gene]['S']].remove(dying_gene)

            # for duplications we edit the extant list for the dying_gene species
            for born_gene in born_genes:
                # if _is_lost(iTree, born_gene):
                    # continue
                iTree.graph['extants'][iTree.node[dying_gene]['S']].add(born_gene)

        else:
            # for speciation, if the species extants sets don't exist, we initialise them
            if iTree.node[born_genes[0]]['S'] not in iTree.graph['extants']:
                for born_gene in born_genes:
                    # if _is_lost(iTree,born_gene):
                        # continue
                    iTree.graph['extants'][iTree.node[born_gene]['S']].update(
                        iTree.graph['extants'][iTree.node[dying_gene]['S']])

            # now both these lists will be updated to reflect the branched gene
            for born_gene in born_genes:

                # if _is_lost(iTree,born_gene):
                    # continue
                try:
                    iTree.graph['extants'][iTree.node[born_gene]['S']].remove(dying_gene)
                except Exception as e:
                    pass
                    # print('born_gene  ', born_gene, iTree.node[born_gene])
                    # print('dying_gene ', dying_gene, iTree.node[dying_gene])
                    # for key in iTree.graph['extants']:
                    #     print(key, iTree.graph['extants'][key])
                    # raise e

                # if _is_lost(iTree,born_gene):
                    # continue
                iTree.graph['extants'][iTree.node[born_gene]['S']].add(born_gene)

            # finally, the branched gene is longer extant in the ancestral species
            try:
                iTree.graph['extants'][iTree.node[dying_gene]['S']].remove(dying_gene)
            except Exception as e:
                print(born_gene, iTree.node[born_gene])
                print(dying_gene, iTree.node[dying_gene])
                for key in iTree.graph['extants']:
                    print(key, iTree.graph['extants'][key])
                raise e

        # print(dying_gene)
        # print(iTree.graph['extants'])
        # print()
        
        # we now have an up-to-date list of genes coexisting with each of the the born genes
        # we can now add interaction nodes for each of these new genes
        for born_gene in born_genes:
            
            # a lost gene does not have an interactions
            # if _is_lost(iTree,born_gene):
                # continue

            for extant_gene in iTree.graph['extants'][iTree.node[born_gene]['S']]:
                
                interaction = _inode_name(extant_gene, born_gene)

                iTree.add_node(interaction, node_type='interaction', S=iTree.node[born_gene]['S'])

                # we find the existing interaction node that is parent to this novel interaction
                # if the new interaction is a self-interaction, so is the parent interaction
                if born_gene == extant_gene:
                    interaction_parent = _inode_name(dying_gene, dying_gene)
                # if extant_gene also one of the newly born ones, parent is a self interaction
                elif extant_gene in born_genes:
                    interaction_parent = _inode_name(dying_gene, dying_gene)
                # otherwise, the parent is a hetero interaction
                else:
                    interaction_parent = _inode_name(dying_gene, extant_gene)
                
                if interaction_parent in iTree.nodes():
                    # print('born_gene  ', born_gene, iTree.node[born_gene])
                    # print('dying_gene ', dying_gene, iTree.node[dying_gene])
                    # print('interaction ', interaction)
                    # print('interaction_parent ', interaction_parent)
                    # return iTree
                    # raise Exception('missing parent!')
                    iTree.add_edge(interaction_parent, interaction)
    
    iTree.remove_nodes_from([n for n in iTree.nodes()
                            if iTree.node[n].get('node_type', '') == 'gene'])
    
    # clean out all nodes that correspond to lost genes
    for node in iTree.nodes():
        # print('_'.join([gTree.node[g]['name'] for g in node.split('_')]).lower())
        if 'lost' in '_'.join([gTree.node[g]['name'] for g in node.split('_')]).lower():
            iTree.remove_node(node)

    # prune all leaf interactions that don't correspond to extant genes
    extants = set([gTree.node[n]['name'] for n in gTree.nodes() if gTree.out_degree(n) == 0])
    # print(extants)

    superfluous = [n for n in iTree.nodes() if iTree.out_degree(n) == 0 and (
                   gTree.node[n.split('_')[0]]['name'] not in extants or
                   gTree.node[n.split('_')[1]]['name'] not in extants)
                   ]
    # print(len(superfluous))

    # while superfluous:
    #
    #     iTree.remove_nodes_from(superfluous)
    #
    #     superfluous = [n for n in iTree.nodes() if iTree.out_degree(n) == 0 and (
    #                    gTree.node[n.split('_')[0]]['name'] not in extants or
    #                    gTree.node[n.split('_')[1]]['name'] not in extants)
    #                    ]

    # print(len(superfluous))
    # for node in [n for n in iTree.nodes() if iTree.out_degree(n)==0]:
    #     if node.split('_')[0] in extants and node.split('_')[1] in extants:
    #         pass
    #     else:
    #         iTree.remove_node(node)
    
    # strip out all nodes that are not descendents of the new root
    # for node in iTree.nodes():
    #     if not nx.has_path(iTree, _inode_name(root, root), node):
    #         iTree.remove_node(node)

    # all_extants = set()
    # for extants in iTree.graph['extants'].values():
    #     all_extants.update(extants)
    #
    # for node in iTree.nodes():
    #     if set(node.split('_')).intersection(all_extants):
    #         iTree.node[node]['extant'] = True

    return iTree
