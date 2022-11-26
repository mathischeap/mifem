# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/25 6:18 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from objects.mpRfT.base.forms.segment.cochain import mpRfT_SgF_Cochain_Base


class mpRfT2_SgF_Cochain(mpRfT_SgF_Cochain_Base):
    """"""

    def __init__(self, t):
        """"""
        super(mpRfT2_SgF_Cochain, self).__init__(t)
        self._freeze_self_()





if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/segment/base/cochain.py
    pass
