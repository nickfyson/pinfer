# -*- coding: utf-8 -*-

__author__ = 'Nick Fyson'
__email__ = 'mail@nickfyson.co.uk'
__version__ = '0.1.0'

# we add the bbn module to the runtime path
import sys
import os
path_comps = [os.path.dirname(os.path.abspath(__file__)), '..', 'bayesian-belief-networks']
sys.path.insert(1, os.path.sep.join(path_comps))
# added at the start of the path to ensure we supercede any existing installation
# not added to the *very* start as that can cause problems...
import bayesian  # NOQA

from .io import load_notung_nhx  # NOQA
from .itree import build_itree  # NOQA
from .visualise import vis_tree  # NOQA
