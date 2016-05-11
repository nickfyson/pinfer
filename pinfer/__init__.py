# -*- coding: utf-8 -*-

__author__ = 'Nick Fyson'
__email__ = 'mail@nickfyson.co.uk'
__version__ = '0.1.0'

import sys
import os

# we add our submodules at the start of the path to supercede any existing installation
# not added to the *very* start as that can cause problems...

# we add the bbn module to the runtime path
path_comps = [os.path.dirname(os.path.abspath(__file__)), '..', 'bayesian-belief-networks']
sys.path.insert(1, os.path.sep.join(path_comps))
import bayesian  # NOQA

# we add the bbn module to the runtime path
path_comps = [os.path.dirname(os.path.abspath(__file__)), '..', 'bayespy']
sys.path.insert(1, os.path.sep.join(path_comps))
# import bayespy  # NOQA

from .io import load_notung_nhx  # NOQA
from .itree import build_itree  # NOQA
from .visualise import vis_tree  # NOQA
