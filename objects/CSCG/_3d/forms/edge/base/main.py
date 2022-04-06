# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""

from objects.CSCG._3d.forms.base import _3dCSCG_FORM_BASE
from objects.CSCG._3d.forms.edge.base.numbering.main import _3dCSCG_Edge_Numbering
from objects.CSCG._3d.forms.edge.base.matrices import _3dCSCG_Edge_Matrices
from objects.CSCG._3d.forms.edge.base.cochain import _3dCSCG_Edge_Cochain
from objects.CSCG._3d.forms.edge.base.error import _3dCSCG_Edge_Error
from objects.CSCG._3d.forms.edge.base.do import _3dCSCG_Edge_DO
from objects.CSCG._3d.forms.edge.base.num import _3dCSCG_EdgeForm_Num

class _3dCSCG_Edge(_3dCSCG_FORM_BASE, ndim=3):
    """
    This is the parent of all 3d standard edge forms.

    :param mesh:
    :param space:
    :param str orientation: 'inner' or 'outer'.
    :param numbering_parameters: The parameters for the numbering. Including scheme name and other parameters.
        When it is a string, we use it as scheme name and it has not other parameters.
    :type numbering_parameters: dict, str
    :param str name:
    """
    def __init__(self, mesh, space, orientation, numbering_parameters, name):
        super().__init__(mesh, space)
        self._NUM_basis_, self._NUM_basis_components_ = \
            getattr(self.space.num_basis, self.__class__.__name__)
        assert orientation in ('inner', 'outer'), " orientation needs to be 'inner' or 'outer'."
        self._orientation_ = orientation
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_edge_form')
        self.standard_properties.name = name
        self._numbering_ = _3dCSCG_Edge_Numbering(self, numbering_parameters)
        self._matrices_ = _3dCSCG_Edge_Matrices(self)
        self._cochain_ = _3dCSCG_Edge_Cochain(self)
        self._error_ = _3dCSCG_Edge_Error(self)
        self._DO_ = _3dCSCG_Edge_DO(self)
        self._num_ = None

        # self._visualize_ = _3dCSCG_Trace_Visualize(self)
        # self._coboundary_ = _3dCSCG_Trace_Coboundary(self)

    def ___PRIVATE_reset_cache___(self):
        """"""
    @property
    def orientation(self):
        return self._orientation_


    @property
    def numbering(self):
        return self._numbering_

    @property
    def matrices(self):
        return self._matrices_

    @property
    def cochain(self):
        return self._cochain_

    @property
    def error(self):
        return self._error_

    @property
    def do(self):
        return self._DO_

    @property
    def num(self):
        if self._num_ is None:
            self._num_ = _3dCSCG_EdgeForm_Num(self)
        return self._num_
