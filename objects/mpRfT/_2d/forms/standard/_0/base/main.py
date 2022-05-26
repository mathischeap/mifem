# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/19 3:53 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.mpRfT._2d.forms.standard.base.main import mpRfT2_StandardFormBase

from objects.mpRfT._2d.forms.standard._0.base.numbering.main import mpRfT2_S0F_Numbering
from objects.mpRfT._2d.forms.standard._0.base.num import mpRfT2_S0F_Num
from objects.mpRfT._2d.forms.standard._0.base.visualize import mpRfT2_S0F_Visualize
from objects.mpRfT._2d.forms.standard._0.base.error import mpRfT2_S0F_Error
from objects.mpRfT._2d.forms.standard._0.base.matrices.main import mpRfT2_S0F_Matrices
from objects.mpRfT._2d.forms.standard._0.base.coboundary import mpRfT2_S0F_Coboundary
from objects.mpRfT._2d.forms.standard._0.base.reconstruct.main import mpRfT2_S0F_Reconstruct
from objects.mpRfT._2d.forms.standard._0.base.discretize.main import mpRfT2_S0F_Discretize
from objects.mpRfT._2d.forms.standard._0.base.migrate import mpRfT2_S0F_Migrate


class mpRfT2_S0F(mpRfT2_StandardFormBase):
    """"""

    def __init__(self, mesh, hybrid, orientation, numbering_parameters, name):
        """"""
        super(mpRfT2_S0F, self).__init__(mesh, hybrid, orientation, name)
        self.standard_properties.___PRIVATE_add_tag___('mpRfT2_standard_0_form')
        self._k_ = 0

        self._numbering_ = mpRfT2_S0F_Numbering(self, numbering_parameters)
        self._num_ = mpRfT2_S0F_Num(self)
        self._error_ = mpRfT2_S0F_Error(self)

        self._reconstruct_ = mpRfT2_S0F_Reconstruct(self)
        self._discretize_ = mpRfT2_S0F_Discretize(self)
        self._migrate_ = mpRfT2_S0F_Migrate(self)

        self._visualize_ = mpRfT2_S0F_Visualize(self)
        self._matrices_ = mpRfT2_S0F_Matrices(self)
        self._coboundary_ = mpRfT2_S0F_Coboundary(self)


    def ___Pr_check_func___(self, func):
        """"""
        assert func.mesh is self.mesh

        if func.__class__.__name__ == 'mpRfT2_Scalar':
            assert func.ftype in ('standard',), \
                f"mpRfT2_S0F FUNC do not accept func mpRfT2_Scalar of ftype {func.ftype}."
        else:
            raise Exception(f"mpRfT2_S0F FUNC do not accept func {func.__class__}")


    @property
    def matrices(self):
        return self._matrices_

    @property
    def visualize(self):
        return self._visualize_

    @property
    def coboundary(self):
        return self._coboundary_



if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/standard/_0/base/main.py
    from __init__ import rfT2

    mesh = rfT2.rm(100)

    mesh.visualize()

