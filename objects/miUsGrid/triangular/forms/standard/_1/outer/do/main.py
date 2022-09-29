# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 7:10 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.miUsGrid.triangular.forms.standard._1.base.do import miUs_Triangular_S1F_Do



class miUs_Triangular_oS1F_Do(miUs_Triangular_S1F_Do):
    """"""

    def __init__(self, sf):
        """"""
        super(miUs_Triangular_oS1F_Do, self).__init__(sf)
        self._freeze_self_()


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
