# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 3:41 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from screws.quadrature import Quadrature
import numpy as np


class miUsTriangular_S2F_Discretize_StandardScalar(FrozenOnly):
    """"""
    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self.___DISCRETIZE_STANDARD_CACHE___ = None
        self._freeze_self_()

    def __call__(self, update_cochain=True, target='CF', quad_degree=None):
        """
        The return cochain is 'locally full local cochain', which means it is mesh-element-wise
        local cochain. So:

        cochainLocal is a dict, whose keys are mesh element numbers, and values (1-d arrays) are
        the local cochains.
        """
        SELF = self._sf_
        nodes = SELF.space.nodes

        p = [SELF.p + 1 for _ in range(2)] if quad_degree is None else quad_degree

        quad_nodes, quad_weights = Quadrature(p).quad
        if self.___DISCRETIZE_STANDARD_CACHE___ is None \
            or quad_degree != self.___DISCRETIZE_STANDARD_CACHE___[0]:
            magic_factor = 0.25
            xi = np.zeros((SELF.num.basis, p[0]+1, p[1]+1))
            et = np.zeros((SELF.num.basis, p[0]+1, p[1]+1))
            volume = np.zeros(SELF.num.basis)
            for j in range(SELF.p):
                for i in range(SELF.p):
                    m = i + j*SELF.p

                    xi[m,...] = (quad_nodes[0][:,np.newaxis].repeat(p[1]+1, axis=1) + 1)\
                              * (nodes[i+1]-nodes[i]
                              )/2 + nodes[i]
                    et[m,...] = (quad_nodes[1][np.newaxis,:].repeat(p[0]+1, axis=0) + 1)\
                              * (nodes[j+1]-nodes[j]
                              )/2 + nodes[j]
                    volume[m] = (nodes[i+1]-nodes[i]) \
                              * (nodes[j+1]-nodes[j])  * magic_factor
            self.___DISCRETIZE_STANDARD_CACHE___ = quad_degree, xi, et, volume
        else:
            xi, et, volume = self.___DISCRETIZE_STANDARD_CACHE___[1:]

        # ----------- target ---------------------------------------------------
        cochainLocal = dict()
        if target == 'CF':
            FUNC = SELF.CF.do.evaluate_func_at_time()[0]
        else:
            raise NotImplementedError(f"I cannot deal with target = {target}.")
        # ============================================================================

        for i in SELF.mesh.elements.indices:
            element = SELF.mesh.elements[i]
            xyz = element.coordinate_transformation.mapping(xi, et)
            detJ = element.coordinate_transformation.Jacobian(xi, et)
            f_xyz = FUNC(*xyz)
            cochainLocal[i] = np.einsum('jkl, k, l, j -> j',
                f_xyz*detJ, quad_weights[0], quad_weights[1],
                volume, optimize='greedy')

        # pass to cochain.local ...
        if update_cochain: SELF.cochain.local = cochainLocal
        # ...
        return 'locally full local cochain', cochainLocal





if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
