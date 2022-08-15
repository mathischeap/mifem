# -*- coding: utf-8 -*-

from screws.freeze.main import FrozenOnly


from objects.CSCG._3d.fields.vector.do.reconstruct.main import _3dCSCG_Vector_Do_Reconstruct
from objects.CSCG._3d.fields.vector.do.cross_product.main import _3dCSCG_Vector_Do_CP
from objects.CSCG._3d.fields.vector.do.inner_product.main import _3dCSCG_Vector_Do_IP

from screws.quadrature import Quadrature


class _3dCSCG_VectorField_DO(FrozenOnly):
    def __init__(self, vf):
        self._vf_ = vf
        self._reconstruct_ = _3dCSCG_Vector_Do_Reconstruct(vf)
        self._cp_ = _3dCSCG_Vector_Do_CP(vf)
        self._ip_ = _3dCSCG_Vector_Do_IP(vf)
        self._freeze_self_()

    def evaluate_func_at_time(self, time=None):
        return self._vf_.___DO_evaluate_func_at_time___(time=time)

    def inner_product_with_space_of(self, other, quad_degree=None):
        return self._vf_.___PRIVATE_do_inner_product_with_space_of___(other, quad_degree=quad_degree)

    @property
    def reconstruct(self):
        return self._reconstruct_

    @property
    def cross_product(self):
        return self._cp_

    @property
    def inner_product(self):
        return self._ip_

    def compute_Ln_norm(self, n=1, quad_degree=None):
        """We compute Ln norm of self.

        int(self)**(n) over the mesh.
        """
        if quad_degree is None:
            quad_degree = (7, 7, 7)

        vf = self._vf_

        if vf.ftype == 'standard':
            quad_nodes, quad_weights = Quadrature(quad_degree).quad_ndim

            print(quad_nodes.shape)


        else:
            raise NotImplementedError(f"")

