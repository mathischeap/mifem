# -*- coding: utf-8 -*-

from components.freeze.main import FrozenOnly






class _3dCSC_Space_Visualize_Matplot(FrozenOnly):
    """"""
    def __init__(self, space):
        self._space_ = space
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        """"""
        raise NotImplementedError()