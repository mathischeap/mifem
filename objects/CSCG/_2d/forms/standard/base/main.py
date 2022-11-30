# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""

from objects.CSCG._2d.forms.base import _2dCSCG_FORM_BASE
from objects.CSCG._2d.forms.standard.base.numbering.main import _2dCSCG_Standard_Form_Numbering
from objects.CSCG._2d.forms.standard.base.do import _2dCSCG_Standard_Form_DO
from objects.CSCG._2d.forms.standard.base.cochain import _2dCSCG_Standard_Form_Cochain
from objects.CSCG._2d.forms.standard.base.coboundary import _2dCSCG_Standard_Form_Coboundary
from objects.CSCG._2d.forms.standard.base.matrices import _2dCSCG_Standard_Form_Matrices
from objects.CSCG._2d.forms.standard.base.error import _2dCSCG_Standard_Form_Error
from objects.CSCG._2d.forms.standard.base.operators.main import _2dCSCG_Standard_Form_Operators
from objects.CSCG._2d.forms.standard.base.dofs.main import _2dCSCG_SF_dofs
from objects.CSCG._2d.forms.standard.base.export.main import _2dCSCG_Standard_Form_Export

from objects.CSCG.base.forms.standard.main import CSCG_Standard_Form
from copy import deepcopy


class _2dCSCG_Standard_Form(CSCG_Standard_Form, _2dCSCG_FORM_BASE, ndim=2):
    """This is the parent of all 2d standard forms.

    :param mesh:
    :param space:
    :param is_hybrid:
    :param orientation:
    :param numbering_parameters: The parameters for the numbering. Including scheme name and other parameters.
        When it is a string, we use it as the scheme name, and it has not other parameters.
    :type numbering_parameters: dict, str
    :param name:
    """
    def __init__(self, mesh, space, is_hybrid, orientation, numbering_parameters, name):
        super(_2dCSCG_Standard_Form, self).__init__(mesh, space, name)
        super(_2dCSCG_Standard_Form, self).___init___()
        self._NUM_basis_, self._NUM_basis_components_ = \
            getattr(self.space.num_basis, self.__class__.__name__)
        assert isinstance(is_hybrid, bool), " isHybrid needs to be bool."
        assert orientation in ('inner', 'outer'), " orientation needs to be 'inner' or 'outer'."
        self._IS_hybrid_ = is_hybrid
        self._orientation_ = orientation
        self.standard_properties.___PRIVATE_add_tag___('2dCSCG_standard_form')
        self.___ARGS___ = (is_hybrid, deepcopy(numbering_parameters), name)

        self._numbering_ = _2dCSCG_Standard_Form_Numbering(self, numbering_parameters)
        self._cochain_ = _2dCSCG_Standard_Form_Cochain(self)
        self._error_ = _2dCSCG_Standard_Form_Error(self)
        self._coboundary_ = _2dCSCG_Standard_Form_Coboundary(self)
        self._matrices_ = _2dCSCG_Standard_Form_Matrices(self)
        self._operators_ = _2dCSCG_Standard_Form_Operators(self)
        self._DO_ = _2dCSCG_Standard_Form_DO(self)
        self._dofs_ = None
        self._export_ = None


    @property
    def shadow(self):
        is_hybrid, numbering_parameters, name = self.___ARGS___
        # noinspection PyArgumentList
        return self.__class__(self.mesh, self.space,
                              is_hybrid=is_hybrid,
                              numbering_parameters=numbering_parameters,
                              name = 'shadow-of-' + name)

    def __repr__(self):
        return f"2dCSCG>{self.orientation}-{self.k}SF>{self.standard_properties.name}:{id(self)}"

    @property
    def numbering(self):
        return self._numbering_

    @property
    def cochain(self):
        return self._cochain_

    @property
    def error(self):
        return self._error_

    @property
    def coboundary(self):
        return self._coboundary_

    @property
    def matrices(self):
        return self._matrices_

    @property
    def operators(self):
        return self._operators_

    @property
    def do(self):
        """If it has too many do methods, we group them in to do."""
        return self._DO_

    @property
    def dofs(self):
        if self._dofs_ is None:
            self._dofs_ = _2dCSCG_SF_dofs(self)
        return self._dofs_

    def ___PRIVATE_do_evaluate_basis_at_meshgrid___(self, xi, eta, compute_xietasigma=True):
        """
        Evaluate the basis functions on ``meshgrid(xi, eta, sigma)``.

        :param xi: A 1d iterable object of floats between -1 and 1.
        :param eta: A 1d iterable object of floats between -1 and 1.
        :type xi: list, tuple, numpy.ndarray
        :type eta: list, tuple, numpy.ndarray
        :param bool compute_xietasigma: (`default`:``True``) If we compute the
            ``meshgrid(xi, eta, sigma, indexing='ij')``.
        :returns: A tuple of outputs:

            1. (None, tuple) -- ``(xi, eta, sigma)`` after ``meshgrid`` and ``ravel('F')``.
            2. tuple -- The evaluated basis functions.
        """
        return self.space.do.evaluate_form_basis_at_meshgrid(
            self.k, xi, eta, compute_xietasigma=compute_xietasigma)

    @property
    def export(self):
        if self._export_ is None:
            self._export_ = _2dCSCG_Standard_Form_Export(self)
        return self._export_