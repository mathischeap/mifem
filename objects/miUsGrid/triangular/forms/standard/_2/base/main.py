# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/20 3:14 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from objects.miUsGrid.triangular.forms.standard.base.main import miUsTriangular_SF_Base

from objects.miUsGrid.triangular.forms.standard._2.base.discretize.main import miUsTriangular_S2F_Discretize
from objects.miUsGrid.triangular.forms.standard._2.base.reconstruct.main import miUsTriangular_S2F_Reconstruct
from objects.miUsGrid.triangular.forms.standard._2.base.do.main import miUs_Triangular_S2F_Do
from objects.miUsGrid.triangular.forms.standard._2.base.matrices.main import miUs_Triangular_S2F_Matrices
from objects.miUsGrid.triangular.forms.standard._2.base.operators.main import miUs_Triangular_S2F_Operators
from objects.miUsGrid.triangular.forms.standard._2.base.error import miUs_Triangular_S2F_Error
from objects.miUsGrid.triangular.forms.standard._2.base.export.main import miUsTriangular_S2F_Export


class miUsTriangular_S2F_Base(miUsTriangular_SF_Base):
    """"""

    def __init__(self, mesh, space, hybrid, orientation, name):
        """"""
        super(miUsTriangular_S2F_Base, self).__init__(mesh, space, hybrid, orientation, 2, name)
        self._discretize_ = miUsTriangular_S2F_Discretize(self)
        self._reconstruct_ = miUsTriangular_S2F_Reconstruct(self)

        self._do_ = miUs_Triangular_S2F_Do(self)
        self._error_ = miUs_Triangular_S2F_Error(self)
        self._matrices_ = miUs_Triangular_S2F_Matrices(self)
        self._operators_ = miUs_Triangular_S2F_Operators(self)
        self._export_ = miUsTriangular_S2F_Export(self)


    def ___Pr_check_CF___(self, func):
        """"""
        assert func.__class__.__name__ == "miUsGrid_Triangular_Scalar", f"I need a miUsGrid_Triangular_Scalar as CF."
        assert func.mesh == self.mesh, f"meshes do not match!"

    @property
    def do(self):
        return self._do_

    @property
    def error(self):
        return self._error_

    @property
    def matrices(self):
        return self._matrices_

    @property
    def operators(self):
        return self._operators_

    @property
    def export(self):
        return self._export_




if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/forms/standard/_2/base/main.py
    from __init__ import miTri

    fc = miTri.call('rand0', 3)

    f2 = fc('2-f-i')





