# -*- coding: utf-8 -*-
"""
Not useful yet!

"""

from components.freeze.main import FrozenOnly


class ProblemBase(FrozenOnly):
    """"""
    def __init__(self):
        self._freeze_self_()

    @property
    def ndim(self):
        raise NotImplementedError()
