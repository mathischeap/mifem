# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/19 3:32 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.main import FrozenClass
from objects.mpRfT.base.forms.base.boundary_condition import mpRfT_Form_BoundaryCondition




class mpRfT_FormBase(FrozenClass):
    """"""

    def __init__(self, mesh, name):
        """"""
        self._mesh_ = mesh
        self._k_ = None
        self.___define_parameters___ = None
        self.standard_properties.___PRIVATE_add_tag___('mpRfT_form')
        self.standard_properties.name = name

        self._analytic_expression_ = None
        self._BC_ = None

    @property
    def mesh(self):
        return self._mesh_

    @property
    def k(self):
        return self._k_

    @property
    def ndim(self):
        """Return the dimensions."""
        return self.mesh.ndim

    @property
    def analytic_expression(self):
        return self._analytic_expression_

    @analytic_expression.setter
    def analytic_expression(self, analytic_expression):
        self.___Pr_check_analytic_expression___(analytic_expression)
        self._analytic_expression_ = analytic_expression

    @property
    def BC(self):
        if self._BC_ is None:
             self._BC_ = mpRfT_Form_BoundaryCondition(self)
        return self._BC_

    #-------- must have methods ------------------------------------------------
    def ___Pr_check_analytic_expression___(self, func):
        raise NotImplementedError()

    def ___Pr_check_BC_analytic_expression___(self, ae):
        raise NotImplementedError()

    @property
    def cochain(self):
        raise NotImplementedError()

    #---------------------------------------------------------------------------------------
    @property
    def ___personal_parameters___(self):
        """
        Personal parameters are parameters that are additional to define_parameters. The full parameters
        are combination of both personal and define parameters.

        IMPORTANT:
        """
        if self.analytic_expression is None:
            AE_parameters = None
        else:
            AE_parameters = self.analytic_expression.___parameters___

        RGW_cochain_local = self.cochain.___Pr_RGW_cochain___()
        # RGW_cochain_local will only be in master, in slaves, it is None.
        return {'AE_current_time': self.analytic_expression.current_time, # to save an object, it must have `current_time`.
                'AE_parameters': AE_parameters,
                'region_wise_cochain_local': RGW_cochain_local}

    @property
    def ___parameters___(self):
        """
        ___parameters___ will be called when you call standard_properties.parameters.

        :return:
        """
        parameters = dict()
        parameters.update(self.___define_parameters___)
        parameters.update(self.___personal_parameters___)
        return parameters







if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
