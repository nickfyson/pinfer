
import networkx as nx
import numpy as np


def annotate_observations(iTree, observations, only_leaves=True):

    # remove all existing observations
    [iTree.node[n].pop('observation', None) for n in iTree.nodes()]

    if only_leaves:
        n_to_annotate = [n for n in iTree.nodes() if not iTree.successors(n)]
    else:
        n_to_annotate = iTree.nodes()

    for node in n_to_annotate:

        observation = observations.get(node, None)

        if not observation:
            continue

        new_node = 'observation_' + node

        CPT = np.zeros((2, 2))
        CPT[0, :] = np.array([1.0, 0.0])
        CPT[1, :] = np.array([0.0, 1.0])

        iTree.add_node(new_node,
                       name=new_node,
                       observation=observation,
                       CPT=CPT,
                       node_type='observation')

        iTree.add_edge(node, new_node)

    # we assume a 50% chance of the ancestral self-interaction existing
    iTree.node[nx.topological_sort(iTree)[0]]['prior'] = [0.5, 0.5]
