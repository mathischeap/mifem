# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/16 2:49 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class nCSCG_RF2_LocalCochainBase(FrozenOnly):
    """Smallest-cell-wise local cochain."""

    def __init__(self, f, LC):
        """"""
        assert '_2nCSCG_RF2_standard_form' in f.standard_properties.tags, f"I need a 2nCSCG RF2 standard form."
        self._f_ = f
        self._mesh_ = f.mesh
        self._LC_ =  self.___Pr_check_and_parse_LC___(LC)
        self._signature_ = self._f_.signature
        self._freeze_self_()

    @property
    def signature(self):
        return self._signature_

    @property
    def ___Pr_is_nCSCG_RF2_Standard_LC___(self):
        return True

    def ___Pr_check_and_parse_LC___(self, LC):
        raise NotImplementedError()



if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
