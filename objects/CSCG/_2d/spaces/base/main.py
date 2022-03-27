# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from screws.freeze.main import FrozenClass
from objects.CSCG._2d.spaces.base.num_basis import NumBasis
from objects.CSCG._2d.spaces.base.local_numbering import LocalNumbering
from objects.CSCG._2d.spaces.base.incidence_matrix import IncidenceMatrix
from objects.CSCG._2d.spaces.base.trace_matrix import TraceMatrix
from objects.CSCG._2d.spaces.base.do import _2dCSCG_space_do








class _2dCSCG_Space(FrozenClass):
    """n-D basis; basis functions."""

    def __init__(self, inputs, ndim):
        if ndim is not None: inputs = [inputs for _ in range(ndim)]
        self._inputs_ = inputs
        self._category_ = None
        self.___PRIVATE_generate_1D_basises___()
        assert self.ndim == 2, " <_2dCSCG_Space> "
        self._num_basis_ = NumBasis(self)
        self._local_numbering_ = LocalNumbering(self)
        self._incidence_matrix_ = IncidenceMatrix(self)
        self._trace_matrix_ = TraceMatrix(self)
        self.___define_parameters___ = None
        self.standard_properties.stamp = '2dCSCG|structured|space'
        self._DO_ = _2dCSCG_space_do(self)
        self._freeze_self_()

    @property
    def ___parameters___(self):
        """ It is mandatory for saving."""
        return self.___define_parameters___

    def __eq__(self, other):
        return self.standard_properties.parameters == other.standard_properties.parameters



    def ___PRIVATE_generate_1D_basises___(self):
        """ """
        ndim = len(self._inputs_)
        basises = ()
        p = ()
        nodes = ()
        for i in range(ndim):
            # noinspection PyUnresolvedReferences
            basises += (self.___1D_basis___(*self._inputs_[i]),)
            p += (basises[i].p,)
            nodes += (basises[i].nodes,)
        self._ndim_, self._basises_, self._p_, self._nodes_ = ndim, basises, p, nodes

    @property
    def do(self):
        return self._DO_

    @property
    def num_basis(self):
        return self._num_basis_

    @property
    def local_numbering(self):
        return self._local_numbering_

    @property
    def incidence_matrix(self):
        return self._incidence_matrix_

    @property
    def trace_matrix(self):
        return self._trace_matrix_




    @property
    def category(self):
        """
        We categorize the FunctionSpace using the particular FunctionSpace name
        `self.__class__.__name__` plus a category called `_category_`.

        `_category_` is from its basises. If basises' categories are the same,
        FunctionSpace will take it. Otherwise, `_category_` is named 'mixed'.

        Returns
        -------
        output : str

        """
        _category_ = []
        for n in range(self.ndim):
            _category_.append(self.basises[n].category)
        _category_ = '|'.join(_category_)
        return self.__class__.__name__ + '-' + _category_

    @property
    def Kronecker(self):
        """
        Return True if the basis functions of this space satisfy the Kronecker
        Delta property. Else, return False.

        Returns
        -------
        output[0] : bool
            If True, then the basis functions of this function space all
            satisfy the Kronecker delta property.
        output[1] : list
            A list of bool, represent if basis functions of `self.basis[i]`
            satisfy the Kronecker delta property.

        """
        _Kronecker_ = [self.basises[n].isKronecker for n in range(self.ndim)]
        return all(_Kronecker_), _Kronecker_

    @property
    def IS_Kronecker(self):
        return self.Kronecker[0]

    @property
    def ndim(self):
        return self._ndim_

    @property
    def basises(self):
        """ Return the 1d basis. """
        return self._basises_

    @property
    def p(self):
        return self._p_

    @property
    def nodes(self):
        return self._nodes_



    def ___PRIVATE_do_evaluate_quadrature___(self, quad_degree, quad_type=None):
        """
        This method is supposed to be over-written in its children. Otherwise, it will
        raise NotImplementedError.
        """
        raise NotImplementedError()