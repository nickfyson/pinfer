# -*- coding: utf-8 -*-

from .initialise import initialise_iTree
from .label import label_birth_death
from .interact import add_all_inodes


def build_itree(gTree):
    """function to construct interaction tree, given suitably annotated gene tree"""

    iTree = initialise_iTree(gTree)

    label_birth_death(iTree)

    add_all_inodes(iTree)

    return iTree


from .interact_original import add_all_inodes as original_add_all_inodes


def original_build_itree(gTree):
    """function to construct interaction tree, given suitably annotated gene tree"""

    iTree = initialise_iTree(gTree)

    label_birth_death(iTree)

    original_add_all_inodes(iTree)

    return iTree
