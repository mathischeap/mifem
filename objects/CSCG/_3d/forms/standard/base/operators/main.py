
from screws.freeze.main import FrozenOnly
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.main import EWC_SparseMatrix
from objects.CSCG._3d.forms.standard.base.operators.helpers.inner import ___Operators_3dCSCG_sf_Inner___



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
        data_generator = ___Operators_3dCSCG_sf_Inner___(self._sf_, other, quad_degree=quad_degree)
        # note that even all mesh elements are unique, we still cache the output because we may use it for multiple times.
        return EWC_SparseMatrix(self._sf_.mesh.elements, data_generator,
                                self._sf_.mesh.elements.___PRIVATE_elementwise_cache_metric_key___)

    def wedge(self, other, quad_degree=None):
        """"""
        data = self._sf_.___PRIVATE_operator_wedge___(other, quad_degree=quad_degree)
        return EWC_SparseMatrix(self._sf_.mesh.elements, data, 'constant')

    def cross_product(self, *args, **kwargs):
        return self._sf_.special.cross_product(*args, **kwargs)