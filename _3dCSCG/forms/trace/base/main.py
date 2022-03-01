# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from _3dCSCG.forms.base import _3dCSCG_FORM_BASE
from _3dCSCG.forms.trace.base.numbering.main import _3dCSCG_Trace_Numbering

from _3dCSCG.forms.trace.base.visualize.main import _3dCSCG_Trace_Visualize

from inheriting.CSCG.form.trace.main_BASE import CSCG_Trace_Form

from _3dCSCG.forms.trace.base.dofs.main import _3dCSCG_Trace_forms_DOFs
from _3dCSCG.forms.trace.base.error import _3dCSCG_Trace_Error
from _3dCSCG.forms.trace.base.matrices import _3dCSCG_Trace_Matrices
from _3dCSCG.forms.trace.base.coboundary import _3dCSCG_Trace_Coboundary
from _3dCSCG.forms.trace.base.cochain import _3dCSCG_Trace_Cochain
from _3dCSCG.forms.trace.base.do import _3dCSCG_Trace_DO


class _3dCSCG_Standard_Trace(CSCG_Trace_Form, _3dCSCG_FORM_BASE, ndim=3):
    """
    This is the parent of all 3d standard trace forms.

    :param mesh:
    :param space:
    :param str orientation: 'inner' or 'outer'.
    :param numbering_parameters: The parameters for the numbering. Including scheme name and other parameters.
        When it is a string, we use it as scheme name, and it has not other parameters.
    :type numbering_parameters: dict, str
    :param str name:
    """
    def __init__(self, mesh, space, orientation, numbering_parameters, name):
        super().__init__(mesh, space)
        self._NUM_basis_, self._NUM_basis_components_, self._NUM_basis_onside_ = \
            getattr(self.space.num_basis, self.__class__.__name__)
        assert orientation in ('inner', 'outer'), " orientation needs to be 'inner' or 'outer'."
        self._orientation_ = orientation
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_trace_form')
        self.standard_properties.name = name
        self._numbering_ = _3dCSCG_Trace_Numbering(self, numbering_parameters)
        self._cochain_ = _3dCSCG_Trace_Cochain(self)

        self._visualize_ = _3dCSCG_Trace_Visualize(self)
        self._matrices_ = _3dCSCG_Trace_Matrices(self)
        self._coboundary_ = _3dCSCG_Trace_Coboundary(self)
        self._dofs_ = None
        self._DO_ = _3dCSCG_Trace_DO(self)
        self._error_ = _3dCSCG_Trace_Error(self)

    def ___PRIVATE_reset_cache___(self):
        self.cochain.___PRIVATE_reset_cache___()
        self.coboundary.___PRIVATE_reset_cache___()

    def ___PRIVATE_do_resemble___(self, obj_or_filename):
        """

        :param obj_or_filename:
        :return:
        """
        raise NotImplementedError()

    @property
    def dofs(self):
        """The dofs of the trace form."""
        if self._dofs_ is None:
            self._dofs_ = _3dCSCG_Trace_forms_DOFs(self)
        return self._dofs_

    @property
    def error(self):
        return self._error_






