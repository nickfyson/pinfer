
import networkx as nx


def analyse_bayespy(tree, update_repeats=20, verbose=False):
    """perform VB inference on networkx representation of interaction tree"""

    from bayespy.nodes import Categorical, Mixture

    variables = {}

    for node in nx.topological_sort(tree):

        # if there are no incoming edges, this node is the initial state
        if tree.in_degree(node) == 0:
            prior = tree.node[node]['prior'].tolist()

            variables[node] = Categorical(prior)

        for s, t in tree.in_edges(node):

            variables[s]

            p_transitions = tree.node[t]['CPT'].tolist()

            variables[target] = Mixture(variables[origin],
                                        Categorical, p_transitions)

    from bayespy.inference import VB

    inf_engine = VB(*variables.values())

    # any node with the property 'observed' has its value fixed
    for node in tree.nodes():
        if 'observed' in tree.node[node]:
            variables[node].observe(tree.node[node]['observed'])

    inf_engine.update(repeat=update_repeats, verbose=verbose)

    for node, variable in variables.items():
        tree.node[node]['belief'] = [variable.pdf(0), variable.pdf(1)]
