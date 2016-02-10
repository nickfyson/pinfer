#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_pinfer
----------------------------------

Tests for `pinfer` module.
"""

import unittest

import networkx as nx
import numpy as np


class TestPolytree(unittest.TestCase):

    def setUp(self):

        tree = nx.DiGraph()

        tree.add_node('R', prior=np.array([0.8, 0.2]))
        tree.add_node('S', prior=np.array([0.9, 0.1]))

        CPT_W = np.zeros((2, 2))
        CPT_W[0, :] = np.array([0.8, 0.2])
        CPT_W[1, :] = np.array([0.0, 1.0])

        tree.add_node('W', CPT=CPT_W)

        CPT_H = np.zeros((2, 2, 2))
        CPT_H[0, 0, :] = np.array([1.0, 0.0])
        CPT_H[0, 1, :] = np.array([0.1, 0.9])
        CPT_H[1, 0, :] = np.array([0.0, 1.0])
        CPT_H[1, 1, :] = np.array([0.0, 1.0])

        tree.add_node('H', CPT=CPT_H)

        tree.add_edge('R', 'W')
        tree.add_edge('R', 'H')
        tree.add_edge('S', 'H')

        self.tree = tree

    def test_polytree(self):

        from pinfer.pinfer import polytree

        polytree(self.tree)

        assert (np.round(self.tree.node['R']['belief'], 8) == np.array([0.8, 0.2])).all()
        assert (np.round(self.tree.node['S']['belief'], 8) == np.array([0.9, 0.1])).all()
        assert (np.round(self.tree.node['W']['belief'], 8) == np.array([0.64, 0.36])).all()
        assert (np.round(self.tree.node['H']['belief'], 8) == np.array([0.728, 0.272])).all()

        self.tree.node['H']['observation'] = np.array([0., 1.])

        polytree(self.tree)

        assert (np.round(self.tree.node['R']['belief'], 8) ==
                np.array([0.26470588, 0.73529412])).all()

        assert (np.round(self.tree.node['S']['belief'], 8) ==
                np.array([0.66176471, 0.33823529])).all()

        assert (np.round(self.tree.node['W']['belief'], 8) ==
                np.array([0.21176471, 0.78823529])).all()

        assert (np.round(self.tree.node['H']['belief'], 8) ==
                np.array([0.00000000, 1.00000000])).all()

        self.tree.node['W']['observation'] = np.array([0., 1.])

        polytree(self.tree)

        assert (np.round(self.tree.node['R']['belief'], 8) ==
                np.array([0.06716418, 0.93283582])).all()

        assert (np.round(self.tree.node['S']['belief'], 8) ==
                np.array([0.83955224, 0.16044776])).all()

        assert (np.round(self.tree.node['W']['belief'], 8) ==
                np.array([0.00000000, 1.00000000])).all()

        assert (np.round(self.tree.node['H']['belief'], 8) ==
                np.array([0.00000000, 1.00000000])).all()

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
