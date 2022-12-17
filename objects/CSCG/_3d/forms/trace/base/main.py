# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from objects.CSCG._3d.forms.base import _3dCSCG_FORM_BASE
from objects.CSCG._3d.forms.trace.base.numbering.main import _3dCSCG_Trace_Numbering

from objects.CSCG.base.forms.trace.main import CSCG_Trace_Form

from objects.CSCG._3d.forms.trace.base.dofs.main import _3dCSCG_Trace_forms_DOFs
from objects.CSCG._3d.forms.trace.base.error import _3dCSCG_Trace_Error
from objects.CSCG._3d.forms.trace.base.matrices import _3dCSCG_Trace_Matrices
from objects.CSCG._3d.forms.trace.base.coboundary import _3dCSCG_Trace_Coboundary
from objects.CSCG._3d.forms.trace.base.cochain import _3dCSCG_Trace_Cochain
from objects.CSCG._3d.forms.trace.base.do import _3dCSCG_Trace_DO


class _3dCSCG_Standard_Trace(CSCG_Trace_Form, _3dCSCG_FORM_BASE, ndim=3):
    """This is the parent of all 3d standard trace forms.

    :param mesh:
    :param space:
    :param str orientation: 'inner' or 'outer'.
    :param numbering_parameters: The parameters for the numbering. Including scheme name and other parameters.
        When it is a string, we use it as scheme name, and it has not other parameters.
    :type numbering_parameters: dict, str
    :param str name:
    """
    def __init__(self, mesh, space, hybrid: bool, orientation: str, numbering_parameters, name):
        super(_3dCSCG_Standard_Trace, self).__init__(mesh, space, name)
        super(_3dCSCG_Standard_Trace, self).___init___()

        _NUM_basis_, self._NUM_basis_components_, self._NUM_basis_onside_ = \
            getattr(self.space.num_basis, self.__class__.__name__)
        self._NUM_basis_ = _NUM_basis_[hybrid]

        assert orientation in ('inner', 'outer'), " orientation needs to be 'inner' or 'outer'."

        self._orientation_ = orientation
        self._hybrid_ = hybrid
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_trace_form')
        self._numbering_ = _3dCSCG_Trace_Numbering(self, numbering_parameters)
        self._cochain_ = None

        self._error_ = None
        self._matrices_ = None
        self._coboundary_ = None
        self._DO_ = None
        self._dofs_ = None

    def __repr__(self):
        return f"3dCSCG>{self.k}TF>{self.standard_properties.name}:{id(self)}"

    @property
    def numbering(self):
        """Collections of all numbering-related sub-properties."""
        return self._numbering_

    @property
    def cochain(self):
        """Collections of all cochain-related sub-properties."""
        if self._cochain_ is None:
            self._cochain_ = _3dCSCG_Trace_Cochain(self)
        return self._cochain_

    @property
    def error(self):
        if self._error_ is None:
            self._error_ = _3dCSCG_Trace_Error(self)
        return self._error_

    @property
    def matrices(self):
        """Collections of all matrices-related sub-properties."""
        if self._matrices_ is None:
            self._matrices_ = _3dCSCG_Trace_Matrices(self)
        return self._matrices_

    @property
    def coboundary(self):
        """Collections of all coboundary-related sub-properties."""
        if self._coboundary_ is None:
            self._coboundary_ = _3dCSCG_Trace_Coboundary(self)
        return self._coboundary_

    @property
    def do(self):
        """Group all methods."""
        if self._DO_ is None:
            self._DO_ = _3dCSCG_Trace_DO(self)
        return self._DO_

    @property
    def dofs(self):
        """ The dofs of the trace form."""
        if self._dofs_ is None:
            self._dofs_ = _3dCSCG_Trace_forms_DOFs(self)
        return self._dofs_
