# -*- coding: utf-8 -*-
from components.freeze.main import FrozenOnly

class _3dCSCG_S3F_DOFs_Matplot(FrozenOnly):
    """"""
    def __init__(self, dofs):
        """"""
        self._dofs_ = dofs
        self._sf_ = dofs._sf_
        self._mesh_ = dofs._sf_.mesh
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        """"""
        raise NotImplementedError()