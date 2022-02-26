# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from root.config import *
from screws.frozen import FrozenOnly
from importlib import import_module
from _2dCSCG.forms.base import _2dCSCG_FORM_BASE
from _2dCSCG.forms.standard.base.numbering.main import _2dCSCG_Standard_Form_Numbering
from _2dCSCG.forms.standard.base.visualize.main import _2dCSCG_FormVisualize

from tools.linear_algebra.elementwise_cache import EWC_SparseMatrix


from inheriting.CSCG.form.standard.main_BASE import CSCG_Standard_Form
from inheriting.CSCG.form.standard.main_BASE import CSCG_Standard_Form_Cochain_BASE
from inheriting.CSCG.form.standard.main_BASE import CSCG_Standard_Form_Coboundary_BASE


class _2dCSCG_Standard_Form(CSCG_Standard_Form, _2dCSCG_FORM_BASE, ndim=2):
    """This is the parent of all 2d standard forms.

    :param mesh:
    :param space:
    :param is_hybrid:
    :param orientation:
    :param numbering_parameters: The parameters for the numbering. Including scheme name and other parameters.
        When it is a string, we use it as scheme name and it has not other parameters.
    :type numbering_parameters: dict, str
    :param name:
    """
    def __init__(self, mesh, space, is_hybrid, orientation, numbering_parameters, name):
        super().__init__(mesh, space)
        self._NUM_basis_, self._NUM_basis_components_ = \
            getattr(self.space.num_basis, self.__class__.__name__)
        assert isinstance(is_hybrid, bool), " isHybrid needs to be bool."
        assert orientation in ('inner', 'outer'), " orientation needs to be 'inner' or 'outer'."
        self._IS_hybrid_ = is_hybrid
        self._orientation_ = orientation
        self.standard_properties.name = name
        self.standard_properties.___PRIVATE_add_tag___('2dCSCG_standard_form')
        self._numbering_ = _2dCSCG_Standard_Form_Numbering(self, numbering_parameters)

        self._cochain_ = _2dCSCG_Standard_Form_Cochain(self)
        self._error_ = _2dCSCG_Standard_Form_Error(self)
        self._coboundary_ = _2dCSCG_Standard_Form_Coboundary(self)
        self._matrices_ = _2dCSCG_Standard_Form_Matrices(self)
        self._operators_ = _2dCSCG_Standard_Form_Operators(self)
        self._visualize_ = _2dCSCG_FormVisualize(self)
        self._DO_ = _2dCSCG_Standard_Form_DO(self)

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
        return self.space.DO_evaluate_form_basis_at_meshgrid(
            self.k, xi, eta, compute_xietasigma=compute_xietasigma)


class _2dCSCG_Standard_Form_DO(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._freeze_self_()

    def RESET_cache(self):
        self._sf_.RESET_cache()

    def evaluate_basis_at_meshgrid(self, *args, **kwargs):
        return self._sf_.___DO_evaluate_basis_at_meshgrid___(*args, **kwargs)

    def discretize(self, *args, **kwargs):
        return self._sf_.discretize(*args, **kwargs)

    def reconstruct(self, *args, **kwargs):
        return self._sf_.reconstruct(*args, **kwargs)



class _2dCSCG_Standard_Form_Cochain(CSCG_Standard_Form_Cochain_BASE):
    def __init__(self, sf):
        super().__init__(sf)

    def local_(self, axis):
        """
        The local cochain along a particular axis.

        :param str axis: The local cochain along which axis? ``x``, ``y``.
        :return: The local cochain dict.
        :rtype: Dict[int, numpy.ndarray]
        """
        assert self._sf_.k not in (0, 2), \
            " <Cochain> : %r is a scalar form; has no idea of axes, use cochain.globis." % self._sf_
        numOfBasisComponents = self._sf_.NUM_basis_components
        localAlongAxis = dict()
        for i in self._sf_.mesh.elements:
            if axis == 'x':
                localAlongAxis[i] = self.local[i][:numOfBasisComponents[0]]
            elif axis == 'y':
                localAlongAxis[i] = self.local[i][numOfBasisComponents[0]:]
            else:
                raise Exception()
        return localAlongAxis



class _2dCSCG_Standard_Form_Error(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._freeze_self_()

    def L(self, n=2, quad_degree=None, upon=False):
        """
        The global :math:`L^2` error; it is global, so slaves first send info to the secretary who computes the
        global error and sends it back to all slaves.

        :param int n: (`default`:``2``) :math:`L^{n}` error.
        :param quad_degree: The quadrature degree used to compute the error.
        :param bool upon: If True, we will shift all the reconstructed value in order to one of the value
            match the exact value. This is very useful to measure the error like the pressure or potential
            which is not determined before-hands. We have to fix it at one point. We can do it before solving,
            but it is very complicated. So we'd better do ``upon=True`` when measure the error.
        :return: The global :math:`L^{n}` error.
        :rtype: float
        """
        assert self._sf_.cochain.local is not None, " I have no cochain."
        OneOrThree = 1 if self._sf_.k in (0, 2) else 2
        quad_degree = [self._sf_.dqp[i] + 2 for i in range(2)] \
            if quad_degree is None else quad_degree
        quad_nodes, _, quad_weights = self._sf_.space.DO_evaluate_quadrature(quad_degree)
        xi, eta = np.meshgrid(*quad_nodes, indexing='ij')
        xi = xi.ravel('F')
        eta = eta.ravel('F')
        xyz, v = self._sf_.reconstruct(*quad_nodes, ravel=True)

        # upon shift ...
        add_to = None
        if upon is True: # we will shift according to the first quadrature point.
            if len(xyz) > 0:
                bi = self._sf_.mesh.elements.indices[0]
                bx, by = xyz[bi]
                bv = v[bi]
                base_xyz = [bx[0], by[0]]
                if OneOrThree == 1:
                    base_val = [bv[0][0],]
                else:
                    base_val = [bv[0][0], bv[1][0]]
            else:
                base_xyz = None
                base_val = None
            base_xyz = cOmm.gather(base_xyz, root=mAster_rank)
            base_val = cOmm.gather(base_val, root=mAster_rank)
            if rAnk == mAster_rank:
                for i, bi in enumerate(base_xyz):
                    if bi is not None:
                        break
                # noinspection PyUnboundLocalVariable
                assert bi is not None
                # noinspection PyUnboundLocalVariable
                base_xyz = bi # find the base_xyz
                # noinspection PyUnboundLocalVariable
                base_val = base_val[i]
                func_val = [self._sf_.func.body[m](*base_xyz) for m in range(OneOrThree)]
                assert len(base_val) == len(func_val)
                add_to = [func_val[m] - base_val[m] for m in range(OneOrThree)]
                # this add_to will be added to all reconstructed value.
            else:
                pass
            add_to = cOmm.bcast(add_to, root=mAster_rank)
        elif upon is False:
            pass
        else:
            raise Exception(f"upon can not be {upon}")
        # ...
        if add_to is not None:
            for i in self._sf_.mesh.elements.indices:
                for j in range(OneOrThree):
                    v[i][j] += add_to[j]
        # ...

        localError = list()
        for i in self._sf_.mesh.elements.indices:
            element = self._sf_.mesh.elements[i]
            detJ = element.coordinate_transformation.Jacobian(xi, eta)
            LEIntermediate = np.sum(
            [(v[i][m] - self._sf_.func.body[m](*xyz[i]))**n for m in range(OneOrThree)], axis=0
            )
            localError.append(np.sum(LEIntermediate * detJ * quad_weights))

        core_local = np.sum(localError)
        core_local = cOmm.gather(core_local, root=mAster_rank)

        if rAnk == mAster_rank:
            globalError = np.sum(core_local) ** (1 / n)
        else:
            globalError = None
        globalError = cOmm.bcast(globalError, root=mAster_rank)

        return globalError

    def H(self, d_func, quad_degree=None, upon=False):
        """
        Global :math:`H^1`-error; :math:`H^1` error includes :math:`H(\mathrm{curl})` and :math:`H(\mathrm{div})`;
        since it basically use the ``globalL2`` method, so the computation is done in the secretary core and
        communications are needed.

        :param d_func: The function of the derivative of ``self`` form.
        :param quad_degree: The quadrature degree used to compute the error.
        :param upon:
        :return: The global :math:`H^{1}` error.
        :rtype: float
        """
        selfErrorL2 = self.L(n=2, quad_degree=quad_degree, upon=upon)
        D_self = self._sf_.coboundary()
        D_self.TW.func.body = d_func
        D_self.TW.current_time = self._sf_.TW.current_time
        D_self.TW.___DO_push_all_to_instant___()
        DErrorL2 = D_self.error.L(n=2, quad_degree=quad_degree)
        return (selfErrorL2 ** 2 + DErrorL2 ** 2) ** 0.5





class _2dCSCG_Standard_Form_Coboundary(CSCG_Standard_Form_Coboundary_BASE):
    def __init__(self, sf):
        super().__init__(sf)

    @property
    def incidence_matrix(self):
        """(scipy.sparse.csc_matrix) Return ths incidence matrix of the standard form."""
        assert self._sf_.k < 2, "volume form has no incidence matrix."
        if self._incidenceMatrix_ is None:
            formName = self._sf_.__class__.__name__
            _incidenceMatrix_ = getattr(self._sf_.space.incidence_matrix, formName)
            self._incidenceMatrix_ = \
                EWC_SparseMatrix(self._sf_.mesh.elements, _incidenceMatrix_, 'constant')
        return self._incidenceMatrix_

    def ___PRIVATE_next_class___(self):
        assert self._sf_.k < 2, "volume form has no next prime space."
        k = self._sf_.k
        nextPath = f'_2dCSCG.forms.standard._{k+1}_form.' + self._sf_.orientation
        nextName = f'_2dCSCG_{k+1}Form_'
        if self._sf_.orientation == 'inner':
            nextName += 'Inner'
        elif self._sf_.orientation == 'outer':
            nextName += 'Outer'
        else:
            raise Exception()
        return getattr(import_module(nextPath), nextName)




class _2dCSCG_Standard_Form_Matrices(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._freeze_self_()

    @property
    def mass(self):
        """(Dict[int, scipy.sparse.csr_matrix]) The mass matrix."""
        return self._sf_.operators.inner(self._sf_)

    @property
    def incidence(self):
        return self._sf_.coboundary.incidence_matrix



class _2dCSCG_Standard_Form_Operators(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._freeze_self_()

    def inner(self, other, quad_degree=None):
        """
        We do ``(self, other)``, and note that here we only return a matrix; we do not do the inner product
        which needs that both forms have cochain.

        :param other: The other form.
        :param quad_degree:
        """
        data_generator = ___Operators_Inner___(self._sf_, other, quad_degree=quad_degree)
        return EWC_SparseMatrix(self._sf_.mesh.elements, data_generator,
                                self._sf_.mesh.elements.___PRIVATE_elementwise_cache_metric_key___)

    def wedge(self, other, quad_degree=None):
        data = self._sf_.___OPERATORS_wedge___(other, quad_degree=quad_degree)
        return EWC_SparseMatrix(self._sf_.mesh.elements, data, 'constant')


class ___Operators_Inner___(FrozenOnly):
    """The class for the inner product matrix."""
    def __init__(self, sf, of, quad_degree=None):
        assert sf.ndim == of.ndim and sf.k == of.k, " <___STORAGE_OPERATORS_INNER___> "
        assert sf.mesh == of.mesh, "Meshes do not match."
        self._mesh_ = sf.mesh
        self._sf_ = sf
        self._of_ = of

        if quad_degree is None:
            quad_degree = [int(np.max([sf.dqp[i], of.dqp[i]])) + 1 for i in range(2)]

        quad_nodes, _, quad_weights = sf.space.DO_evaluate_quadrature(quad_degree, quad_type='Gauss')
        xietasigma, bfSelf = sf.DO.evaluate_basis_at_meshgrid(*quad_nodes)
        if of is sf:
            bfOther = bfSelf
        else:
            xietasigma, bfOther = of.DO.evaluate_basis_at_meshgrid(*quad_nodes)
        self._xietasigma_ = xietasigma
        self._quad_weights_ = quad_weights
        self._bfSelf_ = bfSelf
        self._bfOther_ = bfOther
        self._freeze_self_()

    def __call__(self, i):
        Mi = self._sf_.___OPERATORS_inner___(
            self._of_, i, self._xietasigma_, self._quad_weights_, self._bfSelf_, self._bfOther_
        )
        return Mi