# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/16 12:23 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.nCSCG.rf2.base.form.standard.base.cochain.main import nCSCG_RF2_StandardFormCochainBase


class _2nCSCG_RF2_StandardFormCochain(nCSCG_RF2_StandardFormCochainBase):
    """"""

    def __init__(self, f):
        """"""
        super(_2nCSCG_RF2_StandardFormCochain, self).__init__(f)
        self._freeze_self_()

    def ___Pr_RGW_cochain___(self):
        """"""
        raise NotImplementedError()



if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
