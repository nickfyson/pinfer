#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_pinfer
----------------------------------

Tests for `itree` code.
"""

import unittest


import numpy as np

import networkx as nx

class TestITree(unittest.TestCase):

    def setUp(self):

        from pinfer.io import load_notung_nhx

        self.gTree = load_notung_nhx('tests/data/tree.newick')

        gTree = load_notung_nhx('tests/data/tree.newick')

        from pinfer.itree.evol_time import label_birth_death

        label_birth_death(gTree)

        self.gTree_processed = gTree

    def test_notung_import(self):

        gTree = self.gTree

        assert len(gTree.nodes()) == 383
        assert len(gTree.edges()) == 382

        node = [n for n in gTree.nodes() if gTree.node[n]['name'] == 'n10470'][0]

        dictA = gTree.node[node]
        dictB = {'name': 'n10470', 'S': 'Teleostei', 'D': 'N'}

        assert sorted(dictA.keys()) == sorted(dictB.keys())
        assert sorted(dictA.values()) == sorted(dictB.values())

    def test_gtree_normalisation(self):

        gTree = self.gTree_processed

        # birth time must always be earlier than death for each node...
        for s, t in gTree.edges():
            assert gTree.node[s]['t_death'] == gTree.node[t]['t_birth']

        # birth time of all nodes must be death time of parent...
        for s, t in gTree.edges():
            assert gTree.node[s]['t_death'] == gTree.node[t]['t_birth']

        # all the non-duplication nodes within a given species must be coincident
        for species in set([gTree.node[n]['S'] for n in gTree.nodes()]):
            t_deaths = [np.round(gTree.node[n]['t_death'], 10) for n in gTree.nodes() if (
                        gTree.node[n]['S'] == species and
                        gTree.node[n]['D'] == 'N'
                        )]
        assert len(set(t_deaths)) == 1, 'Not all %s speciations are coincident.' % species

    # the paths from (effective) root to leaf within species must be of unit length
    def test_gtree_branch_lengths(self):

        gTree = self.gTree_processed

        species = {gTree.node[n]['S'] for n in gTree.nodes()}

        trees = {}
        for S in species:

            nodes = [n for n in gTree.nodes() if gTree.node[n]['S']==S]

            parents = []
            for node in nodes:
                parent = gTree.predecessors(node)
                if parent:
                    parents.append(parent[0])

            tree = nx.subgraph(gTree,nodes+parents)

            trees[S] = tree

        for S,tree in trees.items():

            roots  = [n for n in tree.nodes() if not tree.predecessors(n)]
            leaves = [n for n in tree.nodes() if not tree.successors(n)]

            for root in roots:
                for leaf in leaves:

                    for path in nx.all_simple_paths(tree, root, leaf):

                        lengths = [tree.edge[path[i]][path[i+1]]['length'] for i in range(len(path)-1)]

                        residual_length = np.round(sum(lengths),10) - np.round(sum(lengths),0)

                        assert residual_length == 0.0, 'path from %s to %s not of unit length' % \
                                                    (gTree.node[root]['name'], gTree.node[leaf]['name'])


    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
