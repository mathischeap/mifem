# -*- coding: utf-8 -*-



from screws.freeze.base import FrozenOnly
from objects.CSCG._2d.forms.standard.base.operators.helpers.inner import ___Operators_Inner___
from tools.linearAlgebra.elementwiseCache.objects.sparseMatrix.main import EWC_SparseMatrix









class _2dCSCG_Standard_Form_Operators(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._inner_quad_type_ = 'Lobatto'
        self._freeze_self_()

    def inner(self, other, quad_degree=None):
        """
        We do ``(self, other)``, and note that here we only return a matrix; we do not do the inner product
        which needs that both forms have cochain.

        :param other: The other form.
        :param quad_degree:
        """
        data_generator = ___Operators_Inner___(self._sf_, other, quad_degree=quad_degree, quad_type=self._inner_quad_type_)
        return EWC_SparseMatrix(self._sf_.mesh.elements, data_generator,
                                self._sf_.mesh.elements.___PRIVATE_elementwise_cache_metric_key___)

    def wedge(self, other, quad_degree=None):
        data = self._sf_.___PRIVATE_operator_wedge___(other, quad_degree=quad_degree)
        return EWC_SparseMatrix(self._sf_.mesh.elements, data, 'constant')