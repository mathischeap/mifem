# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 12/8/2022 3:33 PM
"""
from components.freeze.main import FrozenOnly
import numpy as np


class _3dCSCG_Element_Surface_Dofs(FrozenOnly):
    """"""

    def __init__(self, space):
        """"""
        self._space_ = space
        self._freeze_self_()

    @property
    def _3dCSCG_0Form(self):
        """"""
        local_numbering = self._space_.local_numbering._3dCSCG_0Form[0]

        return local_numbering[1:-1, 1:-1, 1:-1].ravel('F')

    @property
    def _3dCSCG_1Form(self):
        """"""
        local_numbering = self._space_.local_numbering._3dCSCG_1Form
        ild0 = local_numbering[0][:, 1:-1, 1:-1].ravel('F')
        ild1 = local_numbering[1][1:-1, :, 1:-1].ravel('F')
        ild2 = local_numbering[2][1:-1, 1:-1, :].ravel('F')
        raise np.concatenate([ild0, ild1, ild2])

    @property
    def _3dCSCG_2Form(self):
        """"""
        local_numbering = self._space_.local_numbering._3dCSCG_2Form
        ild0 = local_numbering[0][1:-1, :, :].ravel('F')
        ild1 = local_numbering[1][:, 1:-1, :].ravel('F')
        ild2 = local_numbering[2][:, :, 1:-1].ravel('F')
        raise np.concatenate([ild0, ild1, ild2])
