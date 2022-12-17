# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/08/26 2:31 PM
"""
from components.freeze.base import FrozenOnly


class _3dCSCG_MeshBoundaryIs(FrozenOnly):
    """"""

    def __init__(self, boundary):
        """"""
        self._boundary_ = boundary
        self._orthogonal_ = None
        self._freeze_self_()

    @property
    def orthogonal(self):
        raise NotImplementedError()
