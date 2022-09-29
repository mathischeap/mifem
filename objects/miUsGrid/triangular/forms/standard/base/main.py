# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/20 2:36 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from objects.miUsGrid.triangular.forms.base.main import miUsTriangular_FormBase
from objects.miUsGrid.triangular.forms.standard.base.cochain.main import miUs_Triangular_SF_Cochain
from objects.miUsGrid.triangular.forms.standard.base.num import miUs_Triangular_SF_Num
from objects.miUsGrid.triangular.forms.standard.base.coboundary import miUs_Triangular_SF_Coboundary
from objects.miUsGrid.triangular.forms.standard.base.numbering.main import miUs_Triangular_SF_Numbering

class miUsTriangular_SF_Base(miUsTriangular_FormBase):
    """"""

    def __init__(self, mesh, space, orientation, k, name):
        """"""
        super(miUsTriangular_SF_Base, self).__init__(mesh, space, name)
        assert orientation in ('outer', 'inner'), f"orientation={orientation} invalid."
        self._orientation_ = orientation
        assert k in (0, 1, 2), f"{k}-form is invalid."
        self._k_ = k

        self._cochain_ = miUs_Triangular_SF_Cochain(self)
        self._num_ = miUs_Triangular_SF_Num(self)
        self._coboundary_ = miUs_Triangular_SF_Coboundary(self)
        self._numbering_ = miUs_Triangular_SF_Numbering(self)

    @property
    def orientation(self):
        return self._orientation_

    @property
    def k(self):
        """I am a `k`-form."""
        return self._k_

    @property
    def p(self):
        return self.space.p

    @property
    def coboundary(self):
        return self._coboundary_

    @property
    def numbering(self):
        return self._numbering_


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
