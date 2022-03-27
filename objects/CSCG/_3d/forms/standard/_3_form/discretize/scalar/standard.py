from screws.freeze.base import FrozenOnly
import numpy as np

from screws.quadrature import Quadrature


class _3dCSCG_Discretize_Standard(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self.___DISCRETIZE_STANDARD_CACHE___ = None
        self._freeze_self_()

    def __call__(self, update_cochain:bool=True, quad_degree=None):
        """
        The return cochain is 'locally full local cochain', which means it is mesh-element-wise
        local cochain. So:

        cochainLocal is a dict, whose keys are mesh element numbers, and values (1-d arrays) are
        the local cochains.
        """
        SELF = self._sf_

        p = [SELF.dqp[i] + 1 for i in range(SELF.ndim)] if quad_degree is None else quad_degree
        quad_nodes, quad_weights = Quadrature(p).quad
        if self.___DISCRETIZE_STANDARD_CACHE___ is None \
            or quad_degree != self.___DISCRETIZE_STANDARD_CACHE___[0]:
            magic_factor = 0.125
            xi = np.zeros((SELF.num.basis, p[0]+1, p[1]+1, p[2]+1))
            et = np.zeros((SELF.num.basis, p[0]+1, p[1]+1, p[2]+1))
            si = np.zeros((SELF.num.basis, p[0]+1, p[1]+1, p[2]+1))
            volume = np.zeros(SELF.num.basis)
            for k in range(SELF.p[2]):
                for j in range(SELF.p[1]):
                    for i in range(SELF.p[0]):
                        m = i + j*SELF.p[0] + k*SELF.p[0]*SELF.p[1]
                        xi[m,...] = (quad_nodes[0][:,np.newaxis].repeat(p[1]+1,
                          axis=1)[:,:,np.newaxis].repeat(p[2]+1, axis=2) + 1)\
                                  * (SELF.space.nodes[0][i+1]-SELF.space.nodes[0][i]
                                  )/2 + SELF.space.nodes[0][i]
                        et[m,...] = (quad_nodes[1][np.newaxis,:].repeat(p[0]+1,
                          axis=0)[:,:,np.newaxis].repeat(p[2]+1, axis=2) + 1)\
                                  * (SELF.space.nodes[1][j+1]-SELF.space.nodes[1][j]
                                  )/2 + SELF.space.nodes[1][j]
                        si[m,...] = (quad_nodes[2][np.newaxis,:].repeat(p[1]+1,
                          axis=0)[np.newaxis,:,:].repeat(p[0]+1, axis=0) + 1)\
                                  * (SELF.space.nodes[2][k+1]-SELF.space.nodes[2][k]
                                  )/2 + SELF.space.nodes[2][k]
                        volume[m] = (SELF.space.nodes[0][i+1]-SELF.space.nodes[0][i]) \
                                  * (SELF.space.nodes[1][j+1]-SELF.space.nodes[1][j]) \
                                  * (SELF.space.nodes[2][k+1]-SELF.space.nodes[2][k]) * magic_factor
            self.___DISCRETIZE_STANDARD_CACHE___ = quad_degree, xi, et, si, volume
        else:
            xi, et, si, volume = self.___DISCRETIZE_STANDARD_CACHE___[1:]
        cochainLocal = dict()
        f = SELF.func.body[0]
        JC = dict()
        for i in SELF.mesh.elements.indices:
            element = SELF.mesh.elements[i]
            typeWr2Metric = element.type_wrt_metric.mark
            xyz = element.coordinate_transformation.mapping(xi, et, si)
            if typeWr2Metric in JC:
                detJ = JC[typeWr2Metric]
            else:
                detJ = element.coordinate_transformation.Jacobian(xi, et, si)
                if isinstance(typeWr2Metric, str):
                    JC[typeWr2Metric] = detJ
            fxyz = f(*xyz)
            cochainLocal[i] = np.einsum('jklm, k, l, m, j -> j',
                fxyz*detJ, quad_weights[0], quad_weights[1], quad_weights[2],
                volume, optimize='greedy'
            )
        # isKronecker? ...
        if not SELF.space.IS_Kronecker: raise NotImplementedError()
        # pass to cochain.local ...
        if update_cochain: SELF.cochain.local = cochainLocal
        # ...
        return 'locally full local cochain', cochainLocal
