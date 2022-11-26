# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/19 6:13 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from objects.mpRfT.base.forms.standard.cochain import mpRfT_SF_Cochain_Base


class mpRfT2_SF_Cochain(mpRfT_SF_Cochain_Base):
    """"""

    def __init__(self, f):
        """"""
        super(mpRfT2_SF_Cochain, self).__init__(f)
        self._freeze_self_()




if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/standard/base/cochain.py
    pass
