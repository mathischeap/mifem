# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly
from objects.CSCG._2d.fields.vector.do.reconstruct.main import _2dCSCG_Vector_Do_Reconstruct
from objects.CSCG._2d.fields.vector.do.inner_product.main import _2CSCG_VectorField_InnerProduct

from components.quadrature import Quadrature
from root.config.main import RANK, MASTER_RANK, COMM, np


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

    def compute_Ln_norm(self, n=2, quad_degree=None):
        """We compute Ln norm of self.

        int(self)**(n) over the mesh.
        """
        if quad_degree is None:
            quad_degree = (7, 7)

        vf = self._vf_

        if vf.ftype == 'standard':
            qnx, qny, quad_weights = Quadrature(quad_degree).quad_ndim

            QuadValue = vf.do.reconstruct(qnx, qny)[1]
            detJ = vf.mesh.elements.coordinate_transformation.Jacobian(qnx, qny)

            local_norm = list()
            for i in QuadValue:

                local_norm.append(
                    np.sum(
                        ((QuadValue[i][0])**n+(QuadValue[i][1])**n+(QuadValue[i][2])**n)
                        * quad_weights
                        * detJ[i]
                    )
                )

            local_norm = np.sum(local_norm)

            local_norm = COMM.gather(local_norm, root=MASTER_RANK)
            if RANK == MASTER_RANK:
                local_norm = np.sum(local_norm) ** (1/n)
            else:
                pass

            local_norm = COMM.bcast(local_norm, root=MASTER_RANK)

            return local_norm

        else:
            raise NotImplementedError(f"L^{n}-norm of 2dCSCG scalar is not implemented.")
