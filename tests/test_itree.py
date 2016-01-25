#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_pinfer
----------------------------------

Tests for `itree` code.
"""

import unittest


class TestITree(unittest.TestCase):

    def setUp(self):

        from pinfer.io import load_notung_nhx

        self.gTree = load_notung_nhx('tests/data/tree.newick')

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

        from pinfer.itree.evol_time import add_normalised_edge_lengths, label_birth_death

        gTree = self.gTree

        add_normalised_edge_lengths(gTree)
        label_birth_death(gTree)

        # birth time must always be earlier than death for each node...
        for s, t in gTree.edges():
            assert gTree.node[s]['t_death'] == gTree.node[t]['t_birth']

        # birth time of all nodes must be death time of parent...
        for s, t in gTree.edges():
            assert gTree.node[s]['t_death'] == gTree.node[t]['t_birth']

        # all the non-duplication nodes within a given species must be coincident
        import numpy as np

        for species in set([gTree.node[n]['S'] for n in gTree.nodes()]):
            t_deaths = [np.round(gTree.node[n]['t_death'], 10) for n in gTree.nodes() if (
                        gTree.node[n]['S'] == species and
                        gTree.node[n]['D'] == 'N'
                        )]
        assert len(set(t_deaths)) == 1, 'Not all %s speciations are coincident.' % species

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
