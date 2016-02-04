
import networkx as nx


def annotate_observations(iTree, observations, only_leaves=True):

    # remove all existing observations
    [iTree.node[n].pop('observation', None) for n in iTree.nodes()]

    if only_leaves:
        n_to_annotate = [n for n in iTree.nodes() if not iTree.successors(n)]
    else:
        n_to_annotate = iTree.nodes()

    for node in n_to_annotate:

        observation = observations.get(node, None)

        if observation:
            iTree.node[node]['observation'] = observation

    # we assume a 50% chance of the ancestral self-interaction existing
    iTree.node[nx.topological_sort(iTree)[0]]['prior'] = [0.5, 0.5]
