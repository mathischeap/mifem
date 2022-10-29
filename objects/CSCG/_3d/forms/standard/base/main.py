# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from objects.CSCG._3d.forms.base import _3dCSCG_FORM_BASE
from objects.CSCG._3d.forms.standard.base.numbering.main import _3dCSCG_Standard_Form_Numbering
from objects.CSCG._3d.forms.standard.base.export.main import _3dCSC_Standard_Form_Export
from copy import deepcopy

from objects.CSCG.base.forms.standard.main import CSCG_Standard_Form

from objects.CSCG._3d.forms.standard.base.dofs.main import _3dCSCG_Standard_forms_DOFs
from objects.CSCG._3d.forms.standard.base.operators.main import _3dCSCG_Standard_Form_Operators
from objects.CSCG._3d.forms.standard.base.matrices import _3dCSCG_Standard_Form_Matrices
from objects.CSCG._3d.forms.standard.base.coboundary import _3dCSCG_Standard_Form_Coboundary
from objects.CSCG._3d.forms.standard.base.error import _3dCSCG_Standard_Form_Error
from objects.CSCG._3d.forms.standard.base.cochain.main import _3dCSCG_Standard_Form_Cochain




class _3dCSCG_Standard_Form(CSCG_Standard_Form, _3dCSCG_FORM_BASE, ndim=3):
    """ This is the parent of all 3d standard forms.

    :param mesh:
    :param space:
    :param bool is_hybrid:
    :param str orientation: 'inner' or 'outer'.
    :param numbering_parameters: The parameters for the numbering. Including scheme name and other parameters.
        When it is a string, we use it as the scheme name, and it has not other parameters.
    :type numbering_parameters: dict, str
    :param str name:
    """
    def __init__(self, mesh, space, is_hybrid, orientation, numbering_parameters, name):
        """"""
        super().__init__(mesh, space)
        super().___init___()
        self._NUM_basis_, self._NUM_basis_components_ = \
            getattr(self.space.num_basis, self.__class__.__name__)
        assert isinstance(is_hybrid, bool), " isHybrid needs to be bool."
        assert orientation in ('inner', 'outer'), " orientation needs to be 'inner' or 'outer'."
        self._IS_hybrid_ = is_hybrid
        self._orientation_ = orientation
        self.standard_properties.name = name
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_standard_form')
        self.___ARGS___ = (is_hybrid, orientation, deepcopy(numbering_parameters), name)

        self._numbering_ = _3dCSCG_Standard_Form_Numbering(self, numbering_parameters)
        self._cochain_ = None
        self._error_ = None
        self._coboundary_ = None
        self._matrices_ = None
        self._operators_ = None
        self._DO_ = None

        self._export_ = None
        self._dofs_ = None


    @property
    def shadow(self):
        """A shadow is another me but without cochain."""
        is_hybrid, orientation, numbering_parameters, name = self.___ARGS___
        return self.__class__(self.mesh, self.space,
                              is_hybrid=is_hybrid,
                              orientation=orientation,
                              numbering_parameters=numbering_parameters,
                              name = 'shadow-of-' + name)

    def __repr__(self):
        return f"3dCSCG>{self.k}SF>{self.standard_properties.name}:{id(self)}"

    def RESET_cache(self):
        self.cochain.RESET_cache()
        self.coboundary.RESET_cache()

    @property
    def numbering(self):
        return self._numbering_

    @property
    def cochain(self):
        """"""
        if self._cochain_ is None:
            self._cochain_ = _3dCSCG_Standard_Form_Cochain(self)
        return self._cochain_

    @property
    def error(self):
        """"""
        if self._error_ is None:
            self._error_ = _3dCSCG_Standard_Form_Error(self)
        return self._error_

    @property
    def coboundary(self):
        """"""
        if self._coboundary_ is None:
            self._coboundary_ = _3dCSCG_Standard_Form_Coboundary(self)
        return self._coboundary_

    @property
    def matrices(self):
        if self._matrices_ is None:
            self._matrices_ = _3dCSCG_Standard_Form_Matrices(self)
        return self._matrices_

    @property
    def operators(self):
        if self._operators_ is None:
            self._operators_ = _3dCSCG_Standard_Form_Operators(self)
        return self._operators_

    @property
    def do(self):
        """If it has too many do methods, we group them in to do."""
        return self._DO_

    @property
    def export(self):
        if self._export_ is None:
            self._export_ = _3dCSC_Standard_Form_Export(self)
        return self._export_

    @property
    def dofs(self):
        if self._dofs_ is None:
            self._dofs_ = _3dCSCG_Standard_forms_DOFs(self)
        return self._dofs_





    def ___PRIVATE_do_evaluate_basis_at_meshgrid___(self, xi, eta, sigma, compute_xietasigma=True):
        """Evaluate the basis functions on ``meshgrid(xi, eta, sigma)``.

        :param xi: A 1d iterable object of floats between -1 and 1.
        :param eta: A 1d iterable object of floats between -1 and 1.
        :param sigma: A 1d iterable object of floats between -1 and 1.
        :type xi: list, tuple, numpy.ndarray
        :type eta: list, tuple, numpy.ndarray
        :type sigma: list, tuple, numpy.ndarray
        :param bool compute_xietasigma: (`default`:``True``) If we compute the
            ``meshgrid(xi, eta, sigma, indexing='ij')``.
        :returns: A tuple of outputs:

            1. (None, tuple) -- ``(xi, eta, sigma)`` after ``meshgrid`` and ``ravel('F')``.
            2. tuple -- The evaluated basis functions.
        """
        return self.space.do.evaluate_form_basis_at_meshgrid(
            self.k, xi, eta, sigma, compute_xietasigma=compute_xietasigma)


    def ___PRIVATE_saving_info___(self):
        """"""
        my_info = dict()
        my_info['name'] = self.__class__.__name__
        my_info['mesh_ID'] = str(self.mesh)
        return my_info