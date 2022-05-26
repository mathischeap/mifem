# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/25 6:18 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.mpRfT.base.forms.trace.cochain import mpRfT_TF_Cochain_Base


class mpRfT2_TF_Cochain(mpRfT_TF_Cochain_Base):
    """"""

    def __init__(self, t):
        """"""
        super(mpRfT2_TF_Cochain, self).__init__(t)
        self._freeze_self_()





if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
