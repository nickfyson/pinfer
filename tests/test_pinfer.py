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


def sprinkler_example(analyse_function):

    sprinkler = nx.DiGraph()

    sprinkler.add_node('R', prior=np.array([0.8, 0.2]))

    sprinkler.add_node('S', prior=np.array([0.9, 0.1]))

    sprinkler.add_node('W', CPT=np.zeros((2, 2)))
    sprinkler.add_edge('R', 'W')
    sprinkler.node['W']['CPT'][0, :] = np.array([0.8, 0.2])
    sprinkler.node['W']['CPT'][1, :] = np.array([0.0, 1.0])

    sprinkler.add_node('H', CPT=np.zeros((2, 2, 2)))
    sprinkler.add_edge('R', 'H')
    sprinkler.add_edge('S', 'H')
    sprinkler.node['H']['CPT'][0, 0, :] = np.array([1.0, 0.0])
    sprinkler.node['H']['CPT'][0, 1, :] = np.array([0.1, 0.9])
    sprinkler.node['H']['CPT'][1, 0, :] = np.array([0.0, 1.0])
    sprinkler.node['H']['CPT'][1, 1, :] = np.array([0.0, 1.0])

    analyse_function(sprinkler)

    assert (np.round(sprinkler.node['R']['belief'], 8) == np.array([0.8, 0.2])).all()
    assert (np.round(sprinkler.node['S']['belief'], 8) == np.array([0.9, 0.1])).all()
    assert (np.round(sprinkler.node['W']['belief'], 8) == np.array([0.64, 0.36])).all()
    assert (np.round(sprinkler.node['H']['belief'], 8) == np.array([0.728, 0.272])).all()

    sprinkler.node['H']['observation'] = np.array([0., 1.])

    analyse_function(sprinkler)

    assert (np.round(sprinkler.node['R']['belief'], 8) ==
            np.array([0.26470588, 0.73529412])).all()

    assert (np.round(sprinkler.node['S']['belief'], 8) ==
            np.array([0.66176471, 0.33823529])).all()

    assert (np.round(sprinkler.node['W']['belief'], 8) ==
            np.array([0.21176471, 0.78823529])).all()

    assert (np.round(sprinkler.node['H']['belief'], 8) ==
            np.array([0.00000000, 1.00000000])).all()

    sprinkler.node['W']['observation'] = np.array([0., 1.])

    analyse_function(sprinkler)

    assert (np.round(sprinkler.node['R']['belief'], 8) ==
            np.array([0.06716418, 0.93283582])).all()

    assert (np.round(sprinkler.node['S']['belief'], 8) ==
            np.array([0.83955224, 0.16044776])).all()

    assert (np.round(sprinkler.node['W']['belief'], 8) ==
            np.array([0.00000000, 1.00000000])).all()

    assert (np.round(sprinkler.node['H']['belief'], 8) ==
            np.array([0.00000000, 1.00000000])).all()


def cancer_example(analyse_function):

    #
    # Twardy, C., Nicholson, a, Korb, K., & McNeil, J. (2004).
    # Data mining cardiovascular bayesian networks
    #

    cancer = nx.DiGraph()

    cancer.add_node('M', prior=np.array([0.1, 0.9]))

    cancer.add_node('S', CPT=np.zeros((2, 2)))
    cancer.add_edge('M', 'S')
    cancer.node['S']['CPT'][0, :] = np.array([0.8, 0.2])
    cancer.node['S']['CPT'][1, :] = np.array([0.2, 0.8])

    cancer.add_node('B', CPT=np.zeros((2, 2)))
    cancer.add_edge('M', 'B')
    cancer.node['B']['CPT'][0, :] = np.array([0.95, 0.05])
    cancer.node['B']['CPT'][1, :] = np.array([0.80, 0.20])

    cancer.add_node('C', CPT=np.zeros((2, 2, 2)))
    cancer.add_edge('B', 'C')
    cancer.add_edge('S', 'C')
    cancer.node['C']['CPT'][0, 0, :] = np.array([0.95, 0.05])
    cancer.node['C']['CPT'][0, 1, :] = np.array([0.20, 0.80])
    cancer.node['C']['CPT'][1, 0, :] = np.array([0.20, 0.80])
    cancer.node['C']['CPT'][1, 1, :] = np.array([0.20, 0.80])

    cancer.add_node('H', CPT=np.zeros((2, 2)))
    cancer.add_edge('B', 'H')
    cancer.node['H']['CPT'][0, :] = np.array([0.4, 0.6])
    cancer.node['H']['CPT'][1, :] = np.array([0.2, 0.8])

    cancer.node['C']['observation'] = np.array([0., 1.])

    analyse_function(cancer)

    assert (np.round(cancer.node['M']['belief'], 3) == np.array([0.036, 0.964])).all()
    assert (np.round(cancer.node['S']['belief'], 3) == np.array([0.068, 0.932])).all()
    assert (np.round(cancer.node['B']['belief'], 3) == np.array([0.767, 0.233])).all()
    assert (np.round(cancer.node['C']['belief'], 3) == np.array([0.000, 1.000])).all()
    assert (np.round(cancer.node['H']['belief'], 3) == np.array([0.353, 0.647])).all()


from pinfer.infer import analyse_polytree, analyse_bbn


class TestPolytree(unittest.TestCase):

    def setUp(self):
        pass

    def test_sprinkler(self):
        sprinkler_example(analyse_polytree)

    # def test_cancer(self):
    #     cancer_example(analyse_polytree)

    def tearDown(self):
        pass


class TestBBN(unittest.TestCase):

    def setUp(self):
        pass

    def test_sprinkler(self):
        sprinkler_example(analyse_bbn)

    def test_cancer(self):
        cancer_example(analyse_bbn)

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
