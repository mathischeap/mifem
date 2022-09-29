# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 4:34 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.miUsGrid.triangular.forms.standard.base.do import miUs_Triangular_SF_Do



class miUs_Triangular_S0F_Do(miUs_Triangular_SF_Do):
    """"""

    def __init__(self, sf):
        """"""
        super(miUs_Triangular_S0F_Do, self).__init__(sf)
        self._freeze_self_()


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
