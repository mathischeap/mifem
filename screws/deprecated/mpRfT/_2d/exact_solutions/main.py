# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/17 3:45 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.main import FrozenClass
from objects.mpRfT._2d.exact_solutions.do import mpRfT2_ES_Do


class ExactSolution(FrozenClass):
    """A parent for all exact solution classes on 2dCSCG meshes."""
    def __init__(self, mesh):
        assert mesh.__class__.__name__ == 'mpRfT2_Mesh', "Need a mpRfT2_Mesh."
        self._mesh_ = mesh
        self._status_ = None
        self.___define_parameters___ = None
        self._do_ = mpRfT2_ES_Do(self)
        self._freeze_self_()

    def ___PRIVATE_set_status___(self, status):
        assert status._es_ is self
        self._status_ = status

    @property
    def mesh(self):
        """The mesh."""
        return self._mesh_

    @property
    def status(self):
        """A property that contains all sub-properties/sub-methods about the
        status (variables, parameters, coefficients and so on).
        """
        return self._status_

    @property
    def valid_time(self):
        return self._status_.valid_time

    @property
    def ___parameters___(self):
        return self.___define_parameters___

    def __eq__(self, other):
        return self.standard_properties.parameters == other.standard_properties.parameters

    @property
    def do(self):
        return self._do_


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
