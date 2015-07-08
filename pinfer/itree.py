# -*- coding: utf-8 -*-
"""module docstring here"""

from copy import deepcopy
import networkx as nx


def _inode_name(geneA, geneB):
    # the name of interaction nodes is a concatenation of the two gene names
    # crucially, these are always sorted so the order in which genes are passed is irrelevant
    return '%s_%s' % tuple(sorted((geneA, geneB)))


def get_ordered_nodes(tree):
    """return list of all nodes, ordered by distance from tree
    optionally, also return distances as a dictionary"""
    
    root = nx.topological_sort(tree)[0]

    nodes_dists = []
    for node in tree.nodes():
        nodes_dists.append((nx.shortest_path_length(tree, root, node, weight='weight'), node))

    ordered_nodes = [node for (dist, node) in sorted(nodes_dists)]

    distances = {node: dist for (dist, node) in sorted(nodes_dists)}

    return ordered_nodes, distances


def build_itree(gTree):
    """function to construction interaction tree, given suitably annotated gene tree"""
    # we initialise as a copy of the gTree, so the new tree can be built on the existing structure
    
    iTree = deepcopy(gTree)
    iTree.graph['name'] = 'iTree'
    
    # all existing nodes are annotated as genes
    for node in iTree.nodes():
        iTree.node[node]['node_type'] = 'gene'
        iTree.node[node]['extant']    = False

    # first we build a list of nodes ordered by the distance from the root
    ordered_nodes, distances = get_ordered_nodes(iTree)

    root = ordered_nodes[0]
    # we initialise with the ancestral self-interaction
    interaction = _inode_name(root, root)
    iTree.add_node(interaction, node_type='interaction', S=iTree.node[root]['S'])
    # we begin with a single extant gene
    iTree.node[root]['extant'] = True

    # we now sweep through the original gene tree,
    # moving the effective time forward till we hit each branching point
    # in each iteration we remove the node, and -add interactions from its children
    for parent_gene in ordered_nodes:
        
        # the parent node has branched, and is no longer extant
        iTree.node[parent_gene]['extant'] = False
        
        # the new genes are retrieved from the graph
        new_genes = [n for n in iTree.edge[parent_gene] if iTree.node[n]['node_type'] == 'gene']
        
        # for each in turn, we make the extant and iterate over extant genes to add all potential
        # interactions
        for new_gene in new_genes:
            
            if 'lost' in iTree.node[new_gene]['name'].lower():
                continue
            
            iTree.node[new_gene]['extant'] = True
            
            for gene in [n for n in iTree.nodes() if iTree.node[n].get('extant', False)]:
                
                # genes can only interact if they are in the same species
                if iTree.node[gene]['S'] != iTree.node[new_gene]['S']:
                    continue

                interaction = _inode_name(gene, new_gene)
                iTree.add_node(interaction, node_type='interaction', S=iTree.node[gene]['S'])
                
                # if gene is also one of the new ones
                if gene in new_genes:
                    gene_parent = iTree.predecessors(gene)[0]
                else:
                    gene_parent = gene

                new_gene_parent = iTree.predecessors(new_gene)[0]

                # finally, we need to determine which is the correct parent for
                # the new interaction node if this is a self interaction, the
                # parent is the self-interaction of the parent gene
                if gene == new_gene:
                    interaction_parent = _inode_name(parent_gene, parent_gene)
                else:
                    interaction_parent = _inode_name(gene_parent, new_gene_parent)

                if interaction_parent not in iTree.nodes():
                    print()
                    print('new_gene ', new_gene, '   new_gene_parent ', new_gene_parent)
                    print(iTree.node[new_gene])
                    print(iTree.node[new_gene_parent])
                    print('gene ', gene, '  gene_parent ', gene_parent)
                    print(iTree.node[gene])
                    print(iTree.node[gene_parent])
                    print(interaction)
                    print(interaction_parent)
                    gene_parent = iTree.predecessors(gene)[0]
                    interaction_parent = _inode_name(gene_parent, new_gene_parent)
                    print(gene_parent, iTree.node[gene_parent])
                    print(interaction_parent)
                if interaction_parent not in iTree.nodes():
                    gene_parent = iTree.predecessors(gene_parent)[0]
                    interaction_parent = _inode_name(gene_parent, new_gene_parent)
                    print(gene_parent, iTree.node[gene_parent])
                    print(interaction_parent)
                if interaction_parent not in iTree.nodes():
                    gene_parent = iTree.predecessors(gene_parent)[0]
                    interaction_parent = _inode_name(gene_parent, new_gene_parent)
                    print(gene_parent, iTree.node[gene_parent])
                    print(interaction_parent)
                if interaction_parent not in iTree.nodes():
                    gene_parent = iTree.predecessors(gene_parent)[0]
                    interaction_parent = _inode_name(gene_parent, new_gene_parent)
                    print(gene_parent, iTree.node[gene_parent])
                    print(interaction_parent)
                if interaction_parent not in iTree.nodes():
                    print('halted')
                    raise Exception
                iTree.add_edge(interaction_parent, interaction)
            
    iTree.remove_nodes_from([n for n in iTree.nodes()
                            if iTree.node[n].get('node_type', '') == 'gene'])

    return iTree
