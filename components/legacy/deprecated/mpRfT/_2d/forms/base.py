# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/19 3:31 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from objects.mpRfT.base.forms.base.main import mpRfT_FormBase


class mpRfT2_FormBase(mpRfT_FormBase):
    """"""
    def __init__(self, mesh, name):
        assert mesh.__class__.__name__ == 'mpRfT2_Mesh'
        super(mpRfT2_FormBase, self).__init__(mesh, name)
        self.standard_properties.___PRIVATE_add_tag___('mpRfT2_form')


    #-------- must have methods ------------------------------------------------
    def ___Pr_check_analytic_expression___(self, func):
        raise NotImplementedError()

    def ___Pr_check_BC_analytic_expression___(self, ae):
        raise NotImplementedError()

    @property
    def cochain(self):
        raise NotImplementedError()


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
