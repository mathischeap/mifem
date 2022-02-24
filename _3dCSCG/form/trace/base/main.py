# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from SCREWS.frozen import FrozenOnly
from TOOLS.linear_algebra.elementwise_cache import EWC_SparseMatrix

from _3dCSCG.form.base import _3dCSCG_FORM_BASE
from _3dCSCG.form.trace.base.numbering.main import _3dCSCG_Trace_Numbering

from _3dCSCG.form.trace.base.visualize.main import _3dCSCG_Trace_Visualize

from INHERITING.CSCG.form.trace.main_BASE import CSCG_Trace_Form, CSCG_Trace_Form_Cochain_BASE
from INHERITING.CSCG.form.trace.main_BASE import CSCG_Trace_Form_Coboundary_BASE

from _3dCSCG.form.trace.base.dofs.main import _3dCSCG_Trace_forms_DOFs
from _3dCSCG.form.trace.base.error import _3dCSCG_Trace_Error





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

    def RESET_cache(self):
        self.cochain.RESET_cache()
        self.coboundary.RESET_cache()

    def ___DO_evaluate_basis_at_meshgrid___(self, xi, eta, sigma, compute_xietasigma=True):
        """
        Evaluate the basis functions on ``meshgrid(xi, eta, sigma)``.

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
        return self.space.DO_evaluate_trace_basis_at_meshgrid(
            self.k, xi, eta, sigma, compute_xietasigma=compute_xietasigma)

    def ___DO_resemble___(self, obj_or_filename):
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


class _3dCSCG_Trace_DO(FrozenOnly):
    def __init__(self, tf):
        self._tf_ = tf
        self._freeze_self_()

    def evaluate_basis_at_meshgrid(self, *args, **kwargs):
        return self._tf_.___DO_evaluate_basis_at_meshgrid___(*args, **kwargs)

    def resemble(self, *args, **kwargs):
        return self._tf_.___DO_resemble___(*args, **kwargs)

    def discretize(self, *args, **kwargs):
        return self._tf_.discretize(*args, **kwargs)

    def reconstruct(self, *args, **kwargs):
        return self._tf_.reconstruct(*args, **kwargs)




class _3dCSCG_Trace_Cochain(CSCG_Trace_Form_Cochain_BASE):
    def __init__(self, tf):
        super().__init__(tf)

    def ___local_2_local_TEW___(self):
        """"""
        BO = self._tf_.NUM_basis_onside
        INDICES = [0,]
        for sn in 'NSWEBF':
            # noinspection PyUnresolvedReferences
            INDICES.append(INDICES[-1]+BO[sn])
        _D_ = {'N':0, 'S':1, 'W':2, 'E':3, 'B':4, 'F':5}
        TEW = dict()
        for i in self._tf_.mesh.trace.elements:
            TE = self._tf_.mesh.trace.elements[i]
            CE, CS = TE.CHARACTERISTIC_element, TE.CHARACTERISTIC_side
            i0 = _D_[CS]
            TEW[i] = self.local[CE][INDICES[i0]:INDICES[i0 + 1]]
        self.local_TEW = TEW




class _3dCSCG_Trace_Coboundary(CSCG_Trace_Form_Coboundary_BASE):
    def __init__(self, tf):
        super().__init__(tf)




class _3dCSCG_Trace_Matrices(FrozenOnly):
    def __init__(self, tf):
        self._tf_ = tf
        self._T_ = None
        self._S_ = None
        self._mass_ = None
        self._freeze_self_()

    @property
    def trace(self):
        """Return the trace matrix."""
        if self._T_ is None:
            k = self._tf_.k
            formName = f'_{int(k)}Trace'
            T = getattr(self._tf_.space.trace_matrix, formName)[0]
            self._T_ = \
                EWC_SparseMatrix(self._tf_.mesh.elements, T, 'constant')
        return self._T_

    @property
    def selective(self):
        """Return the selective (mesh-element -> trace element) matrix.

        Like the trace matrix but without minus sign.
        """
        if self._S_ is None:
            k = self._tf_.k
            formName = f'_{int(k)}Trace'
            S = getattr(self._tf_.space.selective_matrix, formName)[0]
            self._S_ = \
                EWC_SparseMatrix(self._tf_.mesh.elements, S, 'constant')
        return self._S_

    @property
    def mass(self):
        """Return the mass matrix. It is a dict, keys are trace-element numbers, and values are the mass
        matrices on the trace elements."""
        if self._mass_ is None:
            self._mass_ = self._tf_.___PRIVATE_generate_TEW_mass_matrices___()
        return self._mass_
