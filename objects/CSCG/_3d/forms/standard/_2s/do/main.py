# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/05 9:55 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from objects.CSCG._3d.forms.standard.base.do import _3dCSCG_Standard_Form_DO


class _3dCSCG_S2F_Do(_3dCSCG_Standard_Form_DO):
    """"""

    def __init__(self, sf):
        """"""
        super(_3dCSCG_S2F_Do, self).__init__(sf)
        self._freeze_self_()




if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
