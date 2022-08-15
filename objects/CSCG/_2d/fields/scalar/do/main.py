# -*- coding: utf-8 -*-
import numpy as np
from screws.quadrature import Quadrature
from screws.freeze.base import FrozenOnly
from root.config.main import cOmm, rAnk, mAster_rank
from objects.CSCG._2d.fields.scalar.do.reconstruct.main import _2dCSCG_Scalr_Do_Reconstruct
from objects.CSCG._2d.fields.scalar.do.cross_product.main import _2dCSCG_SclarField_CrossProduct



class _2dCSCG_ScalarField_DO(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._reconstruct_ = _2dCSCG_Scalr_Do_Reconstruct(sf)
        self._cross_product_ = _2dCSCG_SclarField_CrossProduct(sf)
        self._freeze_self_()

    def evaluate_func_at_time(self, time=None):
        return self._sf_.___DO_evaluate_func_at_time___(time=time)

    @property
    def reconstruct(self):
        return self._reconstruct_

    @property
    def cross_product(self):
        return self._cross_product_

    def compute_Ln_norm(self, n=1, quad_degree=None):
        """We compute Ln norm of self.

        int(self)**(n) over the mesh.
        """
        if quad_degree is None:
            quad_degree = (7, 7)

        sf = self._sf_

        if sf.ftype == 'standard':
            qnx, qny, quad_weights = Quadrature(quad_degree).quad_ndim

            QuadValue = sf.do.reconstruct(qnx, qny)[1]
            detJ = sf.mesh.elements.coordinate_transformation.Jacobian(qnx, qny)

            local_norm = list()
            for i in QuadValue:

                local_norm.append(np.sum((QuadValue[i][0])**n * quad_weights * detJ[i]))

            local_norm = np.sum(local_norm)

            local_norm = cOmm.gather(local_norm, root=mAster_rank)
            if rAnk == mAster_rank:
                local_norm = np.sum(local_norm) ** (1/n)
            else:
                pass

            local_norm = cOmm.bcast(local_norm, root=mAster_rank)

            return local_norm

        else:
            raise NotImplementedError(f"L^{n}-norm of 2dCSCG scalar is not implemented.")