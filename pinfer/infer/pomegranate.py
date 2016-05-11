
import networkx as nx
import numpy as np
from pomegranate import DiscreteDistribution, ConditionalProbabilityTable, State, BayesianNetwork


def analyse_pomegranate(tree, update_repeats=1e5, verbose=False):
    """perform inference using the pomegranate python module"""

    if nx.cycle_basis(tree.to_undirected()):
        raise Exception('Pomegranate can only be used on polytrees! (no cycles when undirected)')

    variables = {}

    for node in nx.topological_sort(tree):

        # if there are no incoming edges, this node is the initial state
        if tree.in_degree(node) == 0:

            prior = tree.node[node]['prior']

            variables[node] = DiscreteDistribution({0: prior[0], 1: prior[1]})

        else:

            parents = sorted(tree.predecessors(node))

            CPT = []
            for indices, value in np.ndenumerate(tree.node[node]['CPT']):
                CPT.append(list(indices) + [value])

            parent_vars = [variables[parent] for parent in parents]

            variables[node] = ConditionalProbabilityTable(CPT, parent_vars)

    states = {}
    for key, variable in variables.items():
        states[key] = State(variable, name=key)

    network = BayesianNetwork()

    network.add_states(states.values())

    for s, t in tree.edges():
        network.add_transition(states[s], states[t])

    network.bake()

    observations = {}
    for node in [n for n in tree.nodes() if 'observation' in tree.node[n]]:
        if (tree.node[node]['observation'] == np.array([0.0, 1.0])).all():
            observations[node] = 1
        elif (tree.node[node]['observation'] == np.array([1.0, 0.0])).all():
            observations[node] = 0
        else:
            raise Exception('Have to have binary observations!')

    beliefs = network.forward_backward(observations, max_iterations=update_repeats)

    for state, belief in zip(network.states, beliefs):

        tree.node[state.name]['belief'] = np.zeros(2)

        tree.node[state.name]['belief'][0] = belief.parameters[0][0]
        tree.node[state.name]['belief'][1] = belief.parameters[0][1]

    return tree
