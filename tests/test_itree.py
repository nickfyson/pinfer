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
        pass

    def test_notung_import(self):

        from pinfer.io import load_notung_nhx

        gTree = load_notung_nhx('tests/data/tree.newick')

        assert len(gTree.nodes()) == 383
        assert len(gTree.edges()) == 382

        node = [n for n in gTree.nodes() if gTree.node[n]['name'] == 'n10470'][0]

        dictA = gTree.node[node]
        dictB = {'name': 'n10470', 'S': 'Teleostei', 'D': 'N'}

        assert sorted(dictA.keys()) == sorted(dictB.keys())
        assert sorted(dictA.values()) == sorted(dictB.values())

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
