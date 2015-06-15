#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_pinfer
----------------------------------

Tests for `pinfer` module.
"""

import unittest

from pinfer import pinfer

import networkx as nx
import numpy as np


class TestPinfer(unittest.TestCase):

    def setUp(self):
        self.tree = nx.DiGraph(name='interaction tree')

        self.tree.add_node(0, prior=0.5)
        self.tree.add_node(1, prior=0.5)
        self.tree.add_node(2, prior=0.5)

        self.tree.add_edge(0, 1,
                           p_transitions=np.array([[0.5, 0.5], [0.1, 0.9]]))
        self.tree.add_edge(0, 2,
                           p_transitions=np.array([[0.5, 0.5], [0.2, 0.8]]))

    def test_full_inference(self):
        self.tree.node[1]['observed'] = 0

        pinfer.analyse_tree(self.tree, verbose=False)

        assert round(self.tree.node[0]['posterior'], 12) == 0.166666666667
        assert round(self.tree.node[1]['posterior'], 12) == 0.75
        assert round(self.tree.node[2]['posterior'], 12) == 0.557506665976

        self.tree.node[1].pop('observed')

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
