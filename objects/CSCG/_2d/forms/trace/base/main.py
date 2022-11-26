# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""

from objects.CSCG._2d.forms.base import _2dCSCG_FORM_BASE
from objects.CSCG._2d.forms.trace.base.numbering.main import _2dCSCG_Trace_Numbering
from objects.CSCG._2d.forms.trace.base.visualize import _2dCSCG_TraceVisualize
from objects.CSCG._2d.forms.trace.base.do import _2dCSCG_Trace_DO
from objects.CSCG._2d.forms.trace.base.cochain import _2dCSCG_Trace_Cochain
from objects.CSCG._2d.forms.trace.base.coboundary import _2dCSCG_TraceCoboundary
from objects.CSCG._2d.forms.trace.base.matrices import _2dCSCG_TraceMatrices

from objects.CSCG.base.forms.trace.main import CSCG_Trace_Form



class _2dCSCG_Standard_Trace(CSCG_Trace_Form, _2dCSCG_FORM_BASE, ndim=2):
    """
    This is the parent of all 2d standard trace forms.

    :param mesh:
    :param space:
    :param bool is_hybrid:
    :param str orientation: 'inner' or 'outer'.
    :param numbering_parameters: The parameters for the numbering. Including scheme name and other parameters.
        When it is a string, we use it as scheme name, and it has not other parameters.
    :type numbering_parameters: dict, str
    :param str name:
    """
    def __init__(self, mesh, space, is_hybrid, orientation, numbering_parameters, name):
        super().__init__(mesh, space)
        super(_2dCSCG_Standard_Trace, self).___init___()

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

    def __repr__(self):
        return f"2dCSCG>{self.k}TF>{self.standard_properties.name}:{id(self)}"

    @property
    def numbering(self):
        """Collections of all numbering-related sub-properties."""
        return self._numbering_

    @property
    def cochain(self):
        """Collections of all cochain-related sub-properties."""
        return self._cochain_

    @property
    def visualize(self):
        """Collections of all visualization-related sub-properties."""
        return self._visualize_

    @property
    def matrices(self):
        """Collections of all matrices-related sub-properties."""
        return self._matrices_

    @property
    def coboundary(self):
        """Collections of all coboundary-related sub-properties."""
        return self._coboundary_

    @property
    def do(self):
        """Group all methods."""
        return self._DO_


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
        return self.space.do.evaluate_trace_basis_at_meshgrid(
            self.k, xi, eta, compute_xietasigma=compute_xietasigma)

