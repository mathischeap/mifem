# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11:22 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.nCSCG.rf2._2d.form.standard._1.base.main import _2nCSCG_RF2_Standard1FormBase


class _2nCSCG_RF2_OuterStandard1Form(_2nCSCG_RF2_Standard1FormBase):
    """"""

    def __init__(self, mesh):
        """"""
        super(_2nCSCG_RF2_OuterStandard1Form, self).__init__(mesh)
        self._freeze_self_()



if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
