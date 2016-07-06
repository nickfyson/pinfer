

def gene_is_lost(iTree, gene):
    """simple function to determine whether the gene has been lost"""
    return 'lost' in iTree.node[gene]['name'].lower()


def get_inode_name(geneA, geneB):
    # the name of interaction nodes is a concatenation of the two gene names
    # crucially, these are always sorted so the order in which genes are passed is irrelevant
    return '%s-%s' % tuple(sorted((geneA, geneB)))

