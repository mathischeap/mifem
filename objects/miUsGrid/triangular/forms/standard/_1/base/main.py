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

    def __init__(self, mesh, space, orientation, name):
        """"""
        super(miUsTriangular_S1F_Base, self).__init__(mesh, space, orientation, 1, name)
        self._error_ = miUs_Triangular_S1F_Error(self)



    def ___Pr_check_CF___(self, func):
        """"""
        assert func.__class__.__name__ == "miUsGrid_Triangular_Vector", f"I need a miUsGrid_Triangular_Vector as CF."
        assert func.mesh == self.mesh, f"meshes do not match!"

    @property
    def error(self):
        return self._error_

if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/forms/standard/_1/base/main.py
    pass

