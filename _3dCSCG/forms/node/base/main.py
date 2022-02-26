# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""

from _3dCSCG.forms.base import _3dCSCG_FORM_BASE


from _3dCSCG.forms.node.base.numbering.main import _3dCSCG_Node_Numbering



class _3dCSCG_Node(_3dCSCG_FORM_BASE, ndim=3):
    """
    This is the parent of all (only 1 although) 3d standard node forms.

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
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_node_form')
        self.standard_properties.name = name
        self._numbering_ = _3dCSCG_Node_Numbering(self, numbering_parameters)

        # self._cochain_ = _3dCSCG_Trace_Cochain(self)
        #
        # self._visualize_ = _3dCSCG_Trace_Visualize(self)
        # self._matrices_ = _3dCSCG_Trace_Matrices(self)
        # self._coboundary_ = _3dCSCG_Trace_Coboundary(self)
        # self._DO_ = _3dCSCG_Trace_DO(self)

    def RESET_cache(self):
        """"""

    @property
    def NUM_basis(self):
        """(int) Return a int which represent the number of basis function one mesh element has."""
        return self._NUM_basis_

    @property
    def NUM_basis_components(self):
        """Return a dict, keys are corner names, values are the num basis (must be 1) on this corner.
        """
        return self._NUM_basis_components_

    @property
    def orientation(self):
        return self._orientation_