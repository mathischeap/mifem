
from screws.freeze.base import FrozenOnly
from objects.CSCG._2d.fields.vector.do.reconstruct.main import _2dCSCG_Vector_Do_Reconstruct
from objects.CSCG._2d.fields.vector.do.inner_product.main import _2CSCG_VectorField_InnerProduct

from screws.quadrature import Quadrature

class _2dCSCG_VectorField_DO(FrozenOnly):
    def __init__(self, vf):
        self._vf_ = vf
        self._reconstruct_ = _2dCSCG_Vector_Do_Reconstruct(vf)
        self._inner_product_ = _2CSCG_VectorField_InnerProduct(vf)
        self._freeze_self_()

    def evaluate_func_at_time(self, time=None):
        return self._vf_.___DO_evaluate_func_at_time___(time=time)

    @property
    def reconstruct(self):
        return self._reconstruct_

    @property
    def inner_product(self):
        return self._inner_product_

    def compute_Ln_norm(self, n=1, quad_degree=None):
        """We compute Ln norm of self.

        int(self)**(n) over the mesh.
        """
        if quad_degree is None:
            quad_degree = (7, 7)

        vf = self._vf_

        if vf.ftype == 'standard':
            quad_nodes, quad_weights = Quadrature(quad_degree).quad_ndim

            print(quad_nodes.shape)


        else:
            raise NotImplementedError(f"")