# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 7:36 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.miUsGrid.triangular.forms.standard.base.error import miUs_Triangular_SF_Error


class miUs_Triangular_S1F_Error(miUs_Triangular_SF_Error):
    """"""

    def __init__(self, sf):
        """"""
        super(miUs_Triangular_S1F_Error, self).__init__(sf)
        self._freeze_self_()


if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/forms/standard/_1/base/error.py
    pass
