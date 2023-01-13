# -*- coding: utf-8 -*-
import numpy as np
from components.quadrature import Quadrature
from components.freeze.base import FrozenOnly
from root.config.main import COMM, RANK, MASTER_RANK
from objects.CSCG._2d.fields.scalar.do.reconstruct.main import _2dCSCG_Scalr_Do_Reconstruct
from objects.CSCG._2d.fields.scalar.do.cross_product.main import _2dCSCG_SclarField_CrossProduct

from objects.CSCG._2d.fields.scalar.helpers.mul import _2dCSCG_ScaMulHelper1


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

    def inner_product(self, other):
        """"""
        if other.__class__.__name__ == '_2dCSCG_ScalarField':
            if self._sf_.ftype == 'standard' and other.ftype == 'standard':
                o_func = other.func[0]
                s_func = self._sf_.func[0]

                x0 = _2dCSCG_ScaMulHelper1(s_func, o_func)

                ip_vector = self._sf_.__class__(
                    self._sf_.mesh,
                    x0,
                    ftype='standard',
                    valid_time=self._sf_.valid_time,
                    name=self._sf_.standard_properties.name + "-inner-product-" + other.standard_properties.name
                )
                return ip_vector
            else:
                raise NotImplementedError()
        else:
            raise NotImplementedError()

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

            local_norm = COMM.gather(local_norm, root=MASTER_RANK)
            if RANK == MASTER_RANK:
                local_norm = np.sum(local_norm) ** (1/n)
            else:
                pass

            local_norm = COMM.bcast(local_norm, root=MASTER_RANK)

            return local_norm

        else:
            raise NotImplementedError(f"L^{n}-norm of 2dCSCG scalar is not implemented.")
