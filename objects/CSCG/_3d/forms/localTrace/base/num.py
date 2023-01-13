# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/27/2022 9:51 AM
"""
import sys

if './' not in sys.path:
    sys.path.append('./')
from components.freeze.main import FrozenOnly


class _3dCSCG_LocalTrace_NUM(FrozenOnly):
    """"""

    def __init__(self, ltf):
        """"""
        self._ltf_ = ltf
        self._basis_ = None
        self._BON_ = None
        self._freeze_self_()

    @property
    def basis(self):
        """"""
        if self._basis_ is None:
            self._basis_ = getattr(
                self._ltf_.space.num_basis, self._ltf_.__class__.__name__
            )[0][self._ltf_.whether.hybrid]
        return self._basis_

    @property
    def basis_onside(self):
        if self._BON_ is None:
            self._BON_ = getattr(
                self._ltf_.space.num_basis, self._ltf_.__class__.__name__
            )[1]
        return self._BON_

    @property
    def dofs(self):
        """How many local dofs in this core."""
        return self._ltf_.numbering.num_local_dofs

    @property
    def global_dofs(self):
        """How many dofs in total in all cores."""
        return self._ltf_.numbering.gathering.global_num_dofs


if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
