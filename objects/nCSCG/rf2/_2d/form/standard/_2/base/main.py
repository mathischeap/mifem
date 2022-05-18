# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11:21 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.nCSCG.rf2._2d.form.standard.base.main import _2nCSCG_RF2_StandardFormBase


class _2nCSCG_RF2_Standard2FormBase(_2nCSCG_RF2_StandardFormBase):
    """"""

    def __init__(self, mesh):
        """"""
        super(_2nCSCG_RF2_Standard2FormBase, self).__init__(mesh)
        self.standard_properties.___PRIVATE_add_tag___('_2nCSCG_RF2_standard_2_form')


    #-------- must have methods ------------------------------------------------
    def ___Pr_check_func___(self, func):
        raise NotImplementedError()




if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
