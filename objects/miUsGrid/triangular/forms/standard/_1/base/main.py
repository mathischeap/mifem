# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/20 3:14 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from objects.miUsGrid.triangular.forms.standard.base.main import miUsTriangular_SF_Base

from objects.miUsGrid.triangular.forms.standard._1.base.error import miUs_Triangular_S1F_Error


class miUsTriangular_S1F_Base(miUsTriangular_SF_Base):
    """"""

    def __init__(self, mesh, space, hybrid, orientation, name):
        """"""
        super(miUsTriangular_S1F_Base, self).__init__(mesh, space, hybrid, orientation, 1, name)
        self._error_ = miUs_Triangular_S1F_Error(self)
        self._do_ = None
        self._matrices_ = None
        self._operators_ = None
        self._IDT_ = None
        self._export_ = None

    @property
    def do(self):
        return self._do_

    @property
    def matrices(self):
        return self._matrices_

    @property
    def operators(self):
        return self._operators_


    def ___Pr_check_CF___(self, func):
        """"""
        assert func.__class__.__name__ == "miUsGrid_Triangular_Vector", f"I need a miUsGrid_Triangular_Vector as CF."
        assert func.mesh == self.mesh, f"meshes do not match!"

    @property
    def error(self):
        return self._error_

    @property
    def IDT(self):
        return self._IDT_

    @property
    def export(self):
        return self._export_







if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/forms/standard/_1/base/main.py
    pass


