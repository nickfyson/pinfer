# -*- coding: utf-8 -*-
"""module docstring here"""

import sys
import networkx as nx
import numpy as np
from copy import deepcopy


def _initialise_polytree(tree):
    # all diagnostic evidence and diagnostic messages are initialised to [1,1]
    for node in tree.nodes():
        tree.node[node]['diagnostic'] = np.array([1.0, 1.0])
        tree.node[node].pop('causal', None)
    for s, t in tree.edges():
        tree.edge[s][t]['diagnostic'] = np.array([1.0, 1.0])
        tree.edge[s][t].pop('causal', None)

    # causal support can now be added using the _update_node function
    for node in nx.topological_sort(tree):
        _update_node(tree, node, debug_message='initialising')

    # to enable lazy inclusion of new evidence
    # we record that the tree is initialised
    tree.graph['initialised'] = True


def _update_node(tree, node, debug_message=''):

    ##########
    # update causal support based on incoming messages
    ##########
    # for ancestor nodes causal support is just the prior probability
    # for all other nodes we recalculate based on messages from parents
    if tree.predecessors(node):
        # matrix multiplication of CPT by each message in turn
        causal = tree.node[node]['CPT']
        # we take the dot project of the causal message with the CPT
        # for each parent in turn (sorted order, as per definition of CPT)
        for parent in sorted(tree.predecessors(node), reverse=False):
            causal = np.tensordot(tree.edge[parent][node]['causal'], causal, axes=[0, 0])
        tree.node[node]['causal'] = causal
    else:
        tree.node[node]['causal'] = np.array(tree.node[node]['prior'])

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

        if (tree.edge[node][child]['causal'] == np.array([0., 0.])).all():
            print('%s is passing non-sensical causal message!' % node)
            print('causal_support = ', causal_support)
            print('diagnostic_support = ', diagnostic_support)
            print('diagnostic_summary = ', diagnostic_summary)
            print(debug_message)
            sys.exit()

    ##########
    # update outgoing diagnostic message
    ##########
    parents = sorted(tree.predecessors(node))
    for i, parent in enumerate(parents):

        # we swap columns in the CPT such that the one corresponding to the
        # targeted parent is in position zero
        # this is because we *don't* want to sum over this column
        CPTcopy = np.swapaxes(deepcopy(tree.node[node]['CPT']), 0, i)
        # we build a list of other parents that does *not* include
        # the target parent, but still in sorted order
        others  = parents[:i] + parents[i + 1:]

        # we now sum over dimension 1 each time which
        # leaves the target parent dimension in tact, as required
        for other in others:
            CPTcopy = np.tensordot(tree.edge[other][node]['causal'], CPTcopy, axes=[0, 1])

        # we now sum over the parent node dimension
        # but we sum over the 2nd dimension, since we are interested
        # in the probability over the values of the *parent*
        diag_message = np.tensordot(tree.node[node]['diagnostic'], CPTcopy, axes=[0, 1])
        # this vector is now the message passed from node to parent
        tree.edge[parent][node]['diagnostic'] = diag_message

    ##########
    # finally, we can now update the belief for this node
    ##########
    tree.node[node]['belief'] = ((tree.node[node]['causal'] * tree.node[node]['diagnostic']) /
                                 sum(tree.node[node]['causal'] * tree.node[node]['diagnostic']))

    # add (temprorary?) code to catch obvious errors that may occur
    if np.isnan(tree.node[node]['belief'][0]) or np.isnan(tree.node[node]['belief'][1]):
        print('\tupdating %s' % node)
        print(tree.node[node]['causal'])
        print(tree.node[node]['diagnostic'])
        print(tree.node[node].get('CPT', None))
        sys.exit()
        pass

    return


def analyse_polytree(tree, pivot_node=None, verbose=False):
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

    ##########
    # if necessary we initialise causal and diagnostic support values in the tree
    ##########
    if not tree.graph.get('initialised', False):
        if verbose:
            print('Initialising...')
        _initialise_polytree(tree)

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
        if verbose:
            print('No new observations found...')
        return tree

    if not pivot_node:
        if verbose:
            print('Finding pivot node...')
            sys.stdout.flush()
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
        if verbose:
            print('...Pivot found!')
            sys.stdout.flush()
    else:
        ordered_nodes = nx.topological_sort(tree)
        change_set = set(tree.nodes())
        change_set.remove(pivot_node)

    ##########
    # first pass - inwards
    # we can now propagate changes through the change_set, toward the pivot node
    ##########
    # the ordering is reversed, and only those in the change_set addressed
    if verbose:
        print('First pass...')
        sys.stdout.flush()
    for i,node in enumerate([n for n in reversed(ordered_nodes) if n in change_set]):
        _update_node(tree, node, debug_message='first pass %05d'%i)

    ##########
    # second pass - outwards
    # we now propagate changes from the pivot node out to all other nodes
    ##########
    if verbose:
        print('Second pass...')
        sys.stdout.flush()
    for node in ordered_nodes:
        _update_node(tree, node, debug_message='second pass')

    # finally, we can strip the 'observation' property from all nodes
    # since this evidence has now been incorporated
    for node in tree.nodes():
        if 'observation' in tree.node[node]:
            tree.node[node].pop('observation')
    if verbose:
        print('Complete!')

    return tree
