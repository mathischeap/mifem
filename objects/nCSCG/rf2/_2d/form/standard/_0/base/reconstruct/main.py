# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/16 6:30 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.nCSCG.rf2._2d.form.standard._0.base.reconstruct.BCW_full import _2nCSCG_RF2_S0F_Reconstruct_BCW_full_LCC


class _2nCSCG_RF2_S0F_Reconstruct(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        """"""
        LCC = self._f_.cochain.local

        if LCC.__class__.__name__ == 'nCSCG_RF2__RCW_Full__LocalCochain':
            rct = _2nCSCG_RF2_S0F_Reconstruct_BCW_full_LCC(self._f_)(*args, **kwargs)

        else:
            raise NotImplementedError()

        return rct

    @property
    def matrices(self):
        """The reconstruction matrices."""
        raise NotImplementedError()



if __name__ == "__main__":
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/form/standard/_0/base/reconstruct/main.py
    pass
