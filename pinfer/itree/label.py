
import networkx as nx


def _bump_zero_distance_children(tree):

    for s, t in [(s, t) for s, t in tree.edges() if tree.edge[s][t]['distance'] == 0.0]:
        tree.edge[s][t]['distance'] = 1e-10


def _label_starter_nodes(tree):

    # we construct a list of all non-leaf species, in tree order
    non_leaf_species = []
    for node in [n for n in nx.topological_sort(tree) if tree.successors(n)]:
        # we are only interested in the speciation nodes (ie. *not* duplication)
        if tree.node[node]['D'] == 'N':
            # need to make sure the list is unique
            if tree.node[node]['S'] not in non_leaf_species:
                non_leaf_species.append(tree.node[node]['S'])

    # the speciations are defined to happen at integer times
    species2time = {species: i + 1.0 for i, species in enumerate(non_leaf_species)}

    # all the leaves are fixed to be at the next availble integer time
    leaf_time = len(species2time) + 1.0

    # every non-duplication node is initialised with a t_death property
    # this includes the leaves
    for node in [n for n in tree.nodes() if tree.node[n]['D'] == 'N']:
        tree.node[node]['t_death'] = species2time.get(tree.node[node]['S'], leaf_time)

    # finally, the root is defined to have a t_death of zero
    tree.node[nx.topological_sort(tree)[0]]['t_death'] = 0.0


def _determine_t_death(tree, target):

    # find the time of parent and the distance from it
    parent = tree.predecessors(target)[0]

    start_dist = tree.edge[parent][target]['distance']
    start_time = tree.node[parent]['t_death']

    # build list of descendants within the same species
    descendants = [n for n in nx.descendants(tree, target)
                   if tree.node[n]['S'] == tree.node[target]['S']]

    # find the most distant descendant with 't_death' label
    distances_times = []
    for node in descendants:
        distance_time = (
            nx.shortest_path_length(tree, source=target, target=node, weight='distance'),
            tree.node[node].get('t_death', None)
        )
        distances_times.append(distance_time)

    # max_dist = max(distances)
    distances_times.sort(key=lambda x: x[0])

    end_dist, end_time = distances_times[-1]

    # t_death for node is between that of parent and descendant
    # proportionate to the distance to each
    t_death = start_time + (end_time - start_time) * (start_dist / (start_dist + end_dist))

    tree.node[target]['t_death'] = t_death


def _add_t_births_and_lengths(tree):

    for s, t in tree.edges():

        tree.node[t]['t_birth'] = tree.node[s]['t_death']

        tree.edge[s][t]['length'] = tree.node[t]['t_death'] - tree.node[t]['t_birth']

    # finally, the root is defined to have a birth time of -1.0
    tree.node[nx.topological_sort(tree)[0]]['t_birth'] = -1.0


def label_birth_death(tree):
    """add birth and death time properties to each node"""

    _bump_zero_distance_children(tree)

    # all the speciciation nodes have pre-defined times
    _label_starter_nodes(tree)

    # now we need to label all remaining nodes
    # best achieved in topological order
    for node in [n for n in nx.topological_sort(tree) if 't_death' not in tree.node[n]]:
        _determine_t_death(tree, node)

    _add_t_births_and_lengths(tree)

    return
