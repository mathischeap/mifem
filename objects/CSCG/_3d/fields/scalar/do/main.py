# -*- coding: utf-8 -*-
import numpy as np

from components.freeze.main import FrozenOnly
from objects.CSCG._3d.fields.scalar.do.reconstruct.main import _3dCSCG_Scalar_Do_Reconstruct
from components.quadrature import Quadrature

from root.config.main import COMM, MASTER_RANK, RANK


class _3dCSCG_ScalarField_DO(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._reconstruct_ = _3dCSCG_Scalar_Do_Reconstruct(sf)
        self._freeze_self_()

    def evaluate_func_at_time(self, time=None):
        return self._sf_.___DO_evaluate_func_at_time___(time=time)

    @property
    def reconstruct(self):
        return self._reconstruct_

    def compute_Ln_norm(self, n=1, quad_degree=None):
        """We compute Ln norm of self.

        if n == 1, we do int_(self)

        """
        if quad_degree is None:
            quad_degree = (5, 5, 5)

        sf = self._sf_

        if sf.ftype == 'standard':
            quad_nodes_weights = Quadrature(quad_degree).quad_ndim

            q_xi, q_et, q_sg = quad_nodes_weights[:3]
            q_weights = quad_nodes_weights[-1]

            R = sf.reconstruct(q_xi, q_et, q_sg)[1]

            detJ = sf.mesh.elements.coordinate_transformation.Jacobian(q_xi, q_et, q_sg)

            element_wise_norm = dict()
            NORM = 0
            for i in sf.mesh.elements:
                norm = np.einsum('ijk, ijk, ijk ->', R[i][0]**n, q_weights, detJ[i], optimize='optimal')

                element_wise_norm[i] = norm

                NORM += norm

            NORM = COMM.gather(NORM, root=MASTER_RANK)

            if RANK == MASTER_RANK:
                NORM = np.sum(NORM)**(1/n)

            NORM = COMM.bcast(NORM, root=MASTER_RANK)

            return NORM

        else:
            raise NotImplementedError(f"")
