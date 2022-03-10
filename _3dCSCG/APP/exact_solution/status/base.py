# -*- coding: utf-8 -*-


from root.config.main import *
from screws.quadrature import Quadrature
from screws.freeze.main import FrozenOnly



class Base(FrozenOnly):
    """
    A base (parent) for all exact solution classes.
    """
    def __init__(self, es):
        self._es_ = es
        self._mesh_ = es.mesh

    # ... valid time: override it to make the exact solution only valid at certain time. ...
    @property
    def valid_time(self):
        return None

    @property
    def mesh(self):
        """
        all exact solutions have mesh, no matter if mesh is useful (general exact solutions not use mesh).

        :return: A ``_3dCSCG_Mesh`` object.
        """
        return self._mesh_

    def ___PRIVATE_compute_L_norm_of___(self, what, time=0, n=2, quad_degree=None):
        """
        We compute the :math:`L^n`-norm of the attribute name ``what``.

        :param str what:
        :param float time:
        :param int n: We compute :math:`L^n`-norm.
        :param quad_degree:
        :return:
        """
        assert self._mesh_ is not None, " <MS> : to compute L2_norm, I need a mesh."
        if quad_degree is None:
            quad_degree = int(np.ceil((500000/self._mesh_.elements.GLOBAL_num)**(1/3)))
            if quad_degree > 20: quad_degree = 20
            quad_degree = (quad_degree, quad_degree, quad_degree)
        _Quadrature_ = Quadrature(quad_degree, category='Gauss')
        quad_nodes, quad_weights = _Quadrature_.quad
        _, v = getattr(self, what).reconstruct(*quad_nodes, time=time)
        quad_nodes = np.meshgrid(*quad_nodes, indexing='ij')
        detJ = self.mesh.elements.coordinate_transformation.Jacobian(*quad_nodes)
        LOCAL = 0
        for i in v:
            local = np.sum([vij**n for vij in v[i]], axis=0)
            LOCAL += np.einsum('ijk, ijk, i, j, k -> ', local, detJ[i], *quad_weights, optimize='optimal')
        GLOBAL = cOmm.allreduce(LOCAL, op=MPI.SUM) ** (1/n)
        return GLOBAL