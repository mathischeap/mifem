# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from SCREWS.frozen import FrozenOnly
from _2dCSCG.form.main import _2dCSCG_FORM_BASE
from _2dCSCG.form.trace.numbering.main import _2dCSCG_Trace_Numbering
from _2dCSCG.form.trace.visualize import _2dCSCG_TraceVisualize

from INHERITING.CSCG.form.trace.main_BASE import CSCG_Trace_Form, CSCG_Trace_Form_Cochain_BASE
from INHERITING.CSCG.form.trace.main_BASE import CSCG_Trace_Form_Coboundary_BASE



class _2dCSCG_Standard_Trace(CSCG_Trace_Form, _2dCSCG_FORM_BASE, ndim=2):
    """
    This is the parent of all 2d standard trace forms.

    :param mesh:
    :param space:
    :param bool is_hybrid:
    :param str orientation: 'inner' or 'outer'.
    :param numbering_parameters: The parameters for the numbering. Including scheme name and other parameters.
        When it is a string, we use it as scheme name and it has not other parameters.
    :type numbering_parameters: dict, str
    :param str name:
    """
    def __init__(self, mesh, space, is_hybrid, orientation, numbering_parameters, name):
        super().__init__(mesh, space)
        self._NUM_basis_, self._NUM_basis_components_, self._NUM_basis_onside_ = \
            getattr(self.space.num_basis, self.__class__.__name__)
        assert isinstance(is_hybrid, bool), " isHybrid needs to be bool."
        assert is_hybrid, "Trace form can only be hybrid."
        assert orientation in ('inner', 'outer'), " orientation needs to be 'inner' or 'outer'."
        self._IS_hybrid_ = is_hybrid
        self._orientation_ = orientation
        self.standard_properties.___PRIVATE_add_tag___('2dCSCG_trace_form')
        self.standard_properties.name = name
        self._numbering_ = _2dCSCG_Trace_Numbering(self, numbering_parameters)
        self._cochain_ = _2dCSCG_Trace_Cochain(self)

        self._visualize_ = _2dCSCG_TraceVisualize(self)
        self._matrices_ = _2dCSCG_TraceMatrices(self)
        self._coboundary_ = _2dCSCG_TraceCoboundary(self)
        self._DO_ = _2dCSCG_Trace_DO(self)

    def RESET_cache(self):
        self.cochain.RESET_cache()
        self.coboundary.RESET_cache()

    def ___DO_evaluate_basis_at_meshgrid___(self, xi, eta, compute_xietasigma=True):
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
        return self.space.DO_evaluate_trace_basis_at_meshgrid(
            self.k, xi, eta, compute_xietasigma=compute_xietasigma)

    def ___DO_resemble___(self, obj_or_filename):
        """

        :param obj_or_filename:
        :return:
        """
        raise NotImplementedError()






class _2dCSCG_Trace_DO(FrozenOnly):
    def __init__(self, tf):
        self._tf_ = tf
        self._freeze_self_()

    def evaluate_basis_at_meshgrid(self, *args, **kwargs):
        return self._tf_.___DO_evaluate_basis_at_meshgrid___(*args, **kwargs)

    def resemble(self, *args, **kwargs):
        return self._tf_.___DO_resemble___(*args, **kwargs)





class _2dCSCG_Trace_Cochain(CSCG_Trace_Form_Cochain_BASE):
    def __init__(self, tf):
        super().__init__(tf)

    def ___local_2_local_TEW___(self):
        """"""
        BO = self._tf_.NUM_basis_onside
        INDICES = [0,]
        for sn in 'UDLR':
            # noinspection PyUnresolvedReferences
            INDICES.append(INDICES[-1]+BO[sn])
        _D_ = {'U':0, 'D':1, 'L':2, 'R':3}
        TEW = dict()
        for i in self._tf_.mesh.trace.elements:
            TE = self._tf_.mesh.trace.elements[i]
            CE, Ce = TE.CHARACTERISTIC_element, TE.CHARACTERISTIC_edge
            i0 = _D_[Ce]
            TEW[i] = self.local[CE][INDICES[i0]:INDICES[i0 + 1]]
        self.local_TEW = TEW





class _2dCSCG_TraceCoboundary(CSCG_Trace_Form_Coboundary_BASE):
    def __init__(self, tf):
        super().__init__(tf)




class _2dCSCG_TraceMatrices(FrozenOnly):
    def __init__(self, tf):
        self._tf_ = tf
        self._N_ = None
        self._freeze_self_()

    @property
    def trace(self):
        return self._tf_.coboundary.trace_matrix