
from screws.frozen import FrozenOnly
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.main import EWC_SparseMatrix

import numpy as np



class _3dCSCG_Standard_Form_Operators(FrozenOnly):
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
        # note that even all mesh elements are unique, we still cache the output because we may use it for multiple times.
        return EWC_SparseMatrix(self._sf_.mesh.elements, data_generator,
                                self._sf_.mesh.elements.___PRIVATE_elementwise_cache_metric_key___)

    def wedge(self, other, quad_degree=None):
        """"""
        data = self._sf_.___PRIVATE_operator_wedge___(other, quad_degree=quad_degree)
        return EWC_SparseMatrix(self._sf_.mesh.elements, data, 'constant')

    def cross_product(self, *args, **kwargs):
        return self._sf_.special.cross_product(*args, **kwargs)





class ___Operators_Inner___(FrozenOnly):
    """The class for the inner product matrix."""
    def __init__(self, sf, of, quad_degree=None):
        assert sf.ndim == of.ndim and sf.k == of.k, " <___STORAGE_OPERATORS_INNER___> "
        assert sf.mesh == of.mesh, "Meshes do not match."
        self._mesh_ = sf.mesh
        self._sf_ = sf
        self._of_ = of
        if quad_degree is None:
            quad_degree = [int(np.max([sf.dqp[i], of.dqp[i]])) for i in range(3)]
        quad_nodes, _, quad_weights = sf.space.___PRIVATE_do_evaluate_quadrature___(quad_degree)
        xietasigma, bfSelf = sf.do.evaluate_basis_at_meshgrid(*quad_nodes)
        if of is sf:
            bfOther = bfSelf
        else:
            xietasigma, bfOther = of.do.evaluate_basis_at_meshgrid(*quad_nodes)
        self._xietasigma_ = xietasigma
        self._quad_weights_ = quad_weights
        self._bfSelf_ = bfSelf
        self._bfOther_ = bfOther
        self._freeze_self_()

    def __call__(self, i):
        Mi = self._sf_.___PRIVATE_operator_inner___(
            self._of_, i, self._xietasigma_, self._quad_weights_, self._bfSelf_, self._bfOther_
        )
        return Mi
