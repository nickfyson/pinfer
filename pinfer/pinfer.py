# -*- coding: utf-8 -*-
"""module docstring here"""

import networkx as nx
import numpy as np
from copy import deepcopy


def analyse_tree(tree, update_repeats=20, verbose=True):
    """perform VB inference on networkx representation of interaction tree"""

    from bayespy.nodes import Categorical, Mixture, Gate

    variables = {}

    for node in nx.topological_sort(tree):

        # if there are no incoming edges, this node is the initial state
        if tree.in_degree(node) == 0:
            prior = tree.node[node]['prior']

            variables[node] = Categorical([1 - prior, prior])

        for origin, target in tree.in_edges(node):

            p_transitions = tree.edge[origin][target]['p_transitions']

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
        tree.node[node]['posterior'] = variable.pdf(1)


def analyse_pymc(tree, samples=1000, burns=500):
    """docstring for analyse_pymc"""
    
    import numpy as np
    import pymc as pm

    nodes = {}
    trans = {}

    for node in nx.topological_sort(tree):

        # if there are no incoming edges, this node is the initial state
        if tree.in_degree(node) == 0:
            prior = tree.node[node]['prior']
            
            nodes[node] = pm.Categorical(node, [1 - prior, prior])

        for s, t in tree.in_edges(node):

            p_transitions = tree.edge[s][t]['p_transitions']
            
            trans[t] = pm.Index('tran_%s' % t, p_transitions, nodes[s])
            
            if 'observed' in tree.node[t]:
                nodes[t] = pm.Categorical(t, trans[t],
                                          observed=True, value=tree.node[t]['observed'])
            else:
                nodes[t] = pm.Categorical(t, trans[t])
    
    model = pm.MCMC(list(nodes.values()) + list(trans.values()))

    model.sample(10000, 5000)

    for node in tree.nodes():
        tree.node[node]['posterior'] = np.mean(model.trace(node)[:])

    model = None

    return tree


def polytree(tree):
    """
    Use the message passing algorithm from Pearl 1982 to calculate
    exact posterior probabilities
    
    Implementation of algorithm outlined in Poet & Schachter 1991

    NB all non-root nodes must have a CPT defined
        CPT is the conditional probability, with axes ordered according sorting of parents
        [a,g,h] = sorted(parents)
        CPT axis 0->a, 1->g, 2->h
        CPT.shape = (2,2,2,2) ie. len of parents + 1
    """
    from itertools import combinations

    def update_node(tree, node):

        ##########
        # update causal support based on incoming messages
        ##########
        # for ancestor nodes causal support is just the prior probability
        # for all other nodes we recalculate based on messages from parents
        if tree.predecessors(node):
            # matrix multiplication of CPT by each message in turn
            causal = tree.node[node]['CPT']
            # dot acts on penultimate axis, hence reversed sorting
            for parent in sorted(tree.predecessors(node), reverse=True):
                causal = np.dot(tree.edge[parent][node]['causal'], causal)
            tree.node[node]['causal'] = causal

        ##########
        # update diagnostic support based on incoming messages
        ##########
        if tree.successors(node):
            diagnostic = tree.node[node].get('evidence',
                                             np.ones(len(tree.node[node]['diagnostic'])))
            for child in tree.successors(node):
                diagnostic = diagnostic * tree.edge[node][child]['diagnostic']
            tree.node[node]['diagnostic'] = diagnostic

        ##########
        # update outgoing causal message
        ##########
        for child in tree.successors(node):
            causal_support     = tree.node[node]['causal']
            diagnostic_support = tree.node[node].get('evidence',
                                                     np.ones(len(tree.node[node]['causal'])))

            diagnostic_summary = np.ones(len(causal_support))
            for child_other in tree.successors(node):
                if child != child_other:
                    diagnostic_summary = (diagnostic_summary *
                                          tree.edge[node][child_other]['diagnostic'])
            tree.edge[node][child]['causal'] = (causal_support *
                                                diagnostic_support * diagnostic_summary)
            # tree.edge[node][child]['causal'] = np.random.random((2,))

        ##########
        # update outgoing diagnostic message
        ##########
        parents = sorted(deepcopy(tree.predecessors(node)))
        for i, parent in enumerate(parents):

            CPT     = deepcopy(tree.node[node]['CPT'])

            others  = parents[:i] + parents[i + 1:]
            CPTcopy = np.swapaxes(CPT, 0, i)

            for other in reversed(others):
                CPTcopy = np.dot(tree.edge[other][node]['causal'], CPTcopy)

            diag_message = np.dot(tree.node[node]['diagnostic'], CPTcopy.T)

            tree.edge[parent][node]['diagnostic'] = diag_message

        ##########
        # finally, we can now update the belief for this node
        ##########
        tree.node[node]['belief'] = ((tree.node[node]['causal'] * tree.node[node]['diagnostic']) /
                                     sum(tree.node[node]['causal'] * tree.node[node]['diagnostic']))

        return

    ##########
    # if necessary we initialise causal and diagnostic support values in the tree
    ##########
    if not tree.graph.get('initialised', False):

        for node in nx.topological_sort(tree):
            
            # all diagnostic evidence and diagnostic messages are initialised to [1,1]
            tree.node[node]['diagnostic'] = np.array([1.0, 1.0])
            for child in tree.successors(node):
                tree.edge[node][child]['diagnostic'] = np.array([1.0, 1.0])
            
            # for ancestor nodes causal support is just the prior probability
            if not tree.predecessors(node):
                causal = np.array(tree.node[node]['prior'])
            # otherwise, causal support can be calculated by reference to that of parents
            else:
                # matrix multiplication of CPT by each message in turn
                causal = tree.node[node]['CPT']
                # dot acts on penultimate axis, hence reversed sorting
                for parent in sorted(tree.predecessors(node), reverse=True):
                    causal = np.dot(tree.edge[parent][node]['causal'], causal)

            # having calculated the causal support, we can store it as a property of the node
            tree.node[node]['causal'] = causal
            # in initialising, causal messages are simply equal to that of the node
            for child in tree.successors(node):
                tree.edge[node][child]['causal'] = causal

            # finally, we can now calculate the initial belief for each node
            tree.node[node]['belief'] = ((tree.node[node]['causal'] *
                                          tree.node[node]['diagnostic']) /
                                         sum(tree.node[node]['causal'] *
                                             tree.node[node]['diagnostic']))
        tree.graph['initialised'] = True
    
    ##########
    # we now use the 'observation' property to set the diagnostic evidence for all nodes
    ##########
    for node in tree.nodes():
        if 'observation' in tree.node[node]:
            tree.node[node]['evidence']   = np.array(tree.node[node]['observation'])
            tree.node[node]['diagnostic'] = tree.node[node]['evidence']
    
    ##########
    # find appropriate pivot node in the network
    ##########
    # find set of all nodes that have an observation
    changed = [n for n in tree.nodes() if 'observation' in tree.node[n]]
    # find all nodes found in all paths between all pairs of nodes
    if len(changed) == 0:
        change_set    = set()
        ordered_nodes = tree.nodes()
    else:
        change_set = set()
        for a, b in combinations(changed, 2):
            change_set.update(nx.shortest_path(tree.to_undirected(), a, b))
        if not change_set:
            change_set = set(changed)
        # choice of pivot node is largely arbitrary, so we choose the 'most ancestral' node
        pivot_node = [n for n in nx.topological_sort(tree) if n in change_set][0]
        # we don't need to address the pivot on the first pass, so remove from change_set
        change_set.remove(pivot_node)
        # we build an ordered list of all nodes, furthest from pivot_node first
        ordered_nodes = []
        for node in tree.nodes():
            dist = nx.shortest_path_length(tree.to_undirected(), pivot_node, node)
            ordered_nodes.append((dist, node))
        ordered_nodes = [n for d, n in sorted(ordered_nodes)]
    
    ##########
    # first pass - inwards
    # we can now propagate changes through the change_set, toward the pivot node
    ##########
    # the ordering is reversed, and only those in the change_set addressed
    for node in [n for n in reversed(ordered_nodes) if n in change_set]:
        print('inwards', node)
        update_node(tree, node)
    
    ##########
    # second pass - outwards
    # we now propagate changes from the pivot node out to all other nodes
    ##########
    for node in ordered_nodes:
        print('outwards', node)
        update_node(tree, node)

    # finally, we can strip the 'observation' property from all nodes
    # since this evidence has now been incorporated
    for node in tree.nodes():
        if 'observation' in tree.node[node]:
            tree.node[node].pop('observation')

    return tree
