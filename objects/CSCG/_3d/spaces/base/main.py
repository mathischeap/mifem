# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from components.freeze.main import FrozenClass
from objects.CSCG._3d.spaces.base.num_basis import NumBasis
from objects.CSCG._3d.spaces.base.local_numbering import LocalNumbering
from objects.CSCG._3d.spaces.base.incidence_matrix import IncidenceMatrix
from objects.CSCG._3d.spaces.base.trace_matrix import TraceMatrix
from objects.CSCG._3d.spaces.base.selective_matrix import SelectiveMatrix
from objects.CSCG._3d.spaces.base.do import _3dCSCG_space_do

from objects.CSCG._3d.spaces.base.visualize.main import _3dCSC_Space_Visualize


class _3dCSCG_Space_Base(FrozenClass):
    """n-D basis; basis functions."""
    def __init__(self, inputs, ndim):
        """If `ndim` is not None, we will repeat the inputs for `ndim` times.

        `ndim` to be None or 3!
        """
        assert ndim is None or ndim == 3, f"ndim={ndim} must be None or 3."
        if ndim is not None: inputs = [inputs for _ in range(ndim)]
        self._inputs_ = inputs
        self._category_ = None
        self.___PRIVATE_generate_1D_basises___()
        assert self.ndim == 3, " <FunctionSpace> "
        self._num_basis_ = NumBasis(self)
        self._local_numbering_ = LocalNumbering(self)
        self._incidence_matrix_ = IncidenceMatrix(self)
        self._trace_matrix_ = TraceMatrix(self)
        self._selective_matrix_ = SelectiveMatrix(self)
        self.___define_parameters___ = None
        self.standard_properties.stamp = '3dCSCG|structured|space'
        self._visualize_ = None
        self._DO_ = _3dCSCG_space_do(self)
        self._freeze_self_()

    def __repr__(self):
        """"""
        return f"3dCSCG-space={self.category}={self.p}"

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
            # for all 1d space, it must accept following input types.
            if isinstance(self._inputs_[i], (list, tuple)):
                # noinspection PyUnresolvedReferences
                basises += (self.___1D_basis___(*self._inputs_[i]),)
            elif isinstance(self._inputs_[i], (int, float, str)):
                # noinspection PyUnresolvedReferences
                basises += (self.___1D_basis___(self._inputs_[i]),)
            else:
                raise Exception(f"cannot accept inputs={self._inputs_[i]}")

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
    def selective_matrix(self):
        return self._selective_matrix_


    @property
    def category(self):
        """
        We categorize the FunctionSpace using the particular FunctionSpace name
        `self.__class__.__name__` plus a category called `_category_`.

        `_category_` is from its basis. If basis's categories are the same,
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
        This method is supposed to be over-written in its children.
        """
        raise NotImplementedError()

    @property
    def visualize(self):
        """Use to visualize the 3d spaces."""
        if self._visualize_ is None:
            self._visualize_ = _3dCSC_Space_Visualize(self)
        return self._visualize_