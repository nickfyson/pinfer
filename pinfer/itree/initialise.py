
from copy import deepcopy


def initialise_iTree(tree):

    # we initialise as a copy of the gTree, so the new tree can be built on the existing structure
    tree = deepcopy(tree)
    tree.graph['name'] = 'iTree'

    # all existing nodes are annotated as genes
    for node in tree.nodes():
        tree.node[node]['node_type'] = 'gene'

    return tree
