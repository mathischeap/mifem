# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from screws.freeze.main import FrozenClass

from objects.CSCG._2d.exact_solutions.visualize.main import ExactSolution_Visualize
from objects.CSCG._2d.exact_solutions.do import ExactSolution_do


class ExactSolution(FrozenClass):
    """A parent for all exact solution classes on 2dCSCG meshes."""
    def __init__(self, mesh):
        assert mesh.__class__.__name__ == '_2dCSCG_Mesh', "Need a 2dCSCG mesh."
        self._mesh_ = mesh
        self._visualize_ = ExactSolution_Visualize(self)
        self._do_ = ExactSolution_do(self)
        self._status_ = None
        self.___define_parameters___ = None
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
    def boundary_condition(self):
        """The boundary condition."""
        raise NotImplementedError()

    @property
    def visualize(self):
        """A wrapper of visualization methods."""
        return self._visualize_

    @property
    def do(self):
        return self._do_

    @property
    def ___parameters___(self):
        return self.___define_parameters___

    def __eq__(self, other):
        return self.standard_properties.parameters == other.standard_properties.parameters