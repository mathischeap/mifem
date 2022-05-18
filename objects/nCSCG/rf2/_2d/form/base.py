# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.nCSCG.rf2.base.form.base.main import nCSCG_RF2_FormBase


class _2nCSCG_RF2_FormBase(nCSCG_RF2_FormBase):
    """"""
    def __init__(self, mesh, name):
        super(_2nCSCG_RF2_FormBase, self).__init__(mesh, name)
        self.standard_properties.___PRIVATE_add_tag___('_2nCSCG_RF2_form')


    #-------- must have methods ------------------------------------------------
    def ___Pr_check_func___(self, func):
        raise NotImplementedError()

    @property
    def cochain(self):
        raise NotImplementedError()










if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
