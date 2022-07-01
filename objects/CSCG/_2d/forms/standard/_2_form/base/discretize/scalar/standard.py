# -*- coding: utf-8 -*-
from screws.freeze.base import FrozenOnly

from screws.quadrature import Quadrature
import numpy as np




class _2dCSCG_S2F_Discretize_StandardScalar(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self.___DISCRETIZE_STANDARD_CACHE___ = None
        self._freeze_self_()

    def __call__(self, update_cochain:bool=True, target='func', quad_degree=None):
        """
        The return cochain is 'locally full local cochain', which means it is mesh-element-wise
        local cochain. So:

        cochainLocal is a dict, whose keys are mesh element numbers, and values (1-d arrays) are
        the local cochains.

        :param update_cochain:
        :param quad_degree:
        :return:
        """
        SELF = self._sf_

        p = [SELF.dqp[i] + 1 for i in range(SELF.ndim)] if quad_degree is None else quad_degree
        quad_nodes, quad_weights = Quadrature(p).quad
        if self.___DISCRETIZE_STANDARD_CACHE___ is None \
            or quad_degree != self.___DISCRETIZE_STANDARD_CACHE___[0]:
            magic_factor = 0.25
            xi = np.zeros((SELF.num.basis, p[0]+1, p[1]+1))
            et = np.zeros((SELF.num.basis, p[0]+1, p[1]+1))
            volume = np.zeros(SELF.num.basis)
            for j in range(SELF.p[1]):
                for i in range(SELF.p[0]):
                    m = i + j*SELF.p[0]
                    xi[m,...] = (quad_nodes[0][:,np.newaxis].repeat(p[1]+1, axis=1) + 1)\
                              * (SELF.space.nodes[0][i+1]-SELF.space.nodes[0][i]
                              )/2 + SELF.space.nodes[0][i]
                    et[m,...] = (quad_nodes[1][np.newaxis,:].repeat(p[0]+1, axis=0) + 1)\
                              * (SELF.space.nodes[1][j+1]-SELF.space.nodes[1][j]
                              )/2 + SELF.space.nodes[1][j]
                    volume[m] = (SELF.space.nodes[0][i+1]-SELF.space.nodes[0][i]) \
                              * (SELF.space.nodes[1][j+1]-SELF.space.nodes[1][j])  * magic_factor
            self.___DISCRETIZE_STANDARD_CACHE___ = quad_degree, xi, et, volume
        else:
            xi, et, volume = self.___DISCRETIZE_STANDARD_CACHE___[1:]

        # ----------- target ---------------------------------------------------
        cochainLocal = dict()
        if target == 'func':
            FUNC = SELF.func.body[0]
        else:
            raise NotImplementedError(f"I cannot deal with target = {target}.")
        # ============================================================================

        JC = dict()
        for i in SELF.mesh.elements.indices:
            element = SELF.mesh.elements[i]
            typeWr2Metric = element.type_wrt_metric.mark
            xyz = element.coordinate_transformation.mapping(xi, et)
            if typeWr2Metric in JC:
                detJ = JC[typeWr2Metric]
            else:
                detJ = element.coordinate_transformation.Jacobian(xi, et)
                if isinstance(typeWr2Metric, str):
                    JC[typeWr2Metric] = detJ
            fxyz = FUNC(*xyz)
            cochainLocal[i] = np.einsum('jkl, k, l, j -> j',
                fxyz*detJ, quad_weights[0], quad_weights[1],
                volume, optimize='greedy')

        # isKronecker? ...
        if not SELF.space.IS_Kronecker: raise NotImplementedError()
        # pass to cochain.local ...
        if update_cochain: SELF.cochain.local = cochainLocal
        # ...
        return 'locally full local cochain', cochainLocal