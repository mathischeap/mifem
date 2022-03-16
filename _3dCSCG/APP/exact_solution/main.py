# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from screws.freeze.main import FrozenClass


from _3dCSCG.APP.exact_solution.do import _3dCSCG_ES_DO
from _3dCSCG.APP.exact_solution.visualize import ExactSolution_Visualize


class _3dCSCG_ExactSolution(FrozenClass):
    """A parent for all exact solution classes in _3dCSCG."""
    def __init__(self, mesh):
        assert mesh.__class__.__name__ == '_3dCSCG_Mesh', "Need a 3dCSCG mesh."
        self._mesh_ = mesh
        self._visualize_ = ExactSolution_Visualize(self)
        self._do_ = _3dCSCG_ES_DO(self)
        self._status_ = None
        self.___define_parameters___ = None
        self._freeze_self_()

    def ___PRIVATE_set_status___(self, status):
        assert status._es_ is self
        self._status_ = status

    @property
    def mesh(self):
        """(:class:`_3dCSCG.mesh.mesh._3dCSCG_Mesh`) The mesh."""
        return self._mesh_

    @property
    def status(self):
        """A property that contains all sub-properties/sub-methods about the status(
        variables, parameters, coefficients and so on).
        """
        return self._status_

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