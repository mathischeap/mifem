# -*- coding: utf-8 -*-
import numpy as np
from root.config.main import COMM, MPI
from components.quadrature import Quadrature
from components.freeze.main import FrozenOnly
import random

class Base(FrozenOnly):
    """
    A base (parent) for all exact solution classes.
    """
    def __init__(self, mesh):
        self._mesh_ = mesh

    @property
    def current_time(self):
        """Return a list of current_time of all valid properties.."""
        ct = list()
        for attr_name in self.__dict__:
            attr = getattr(self, attr_name)
            if hasattr(attr, '_current_time_'):
                if attr._current_time_ is not None:
                    ct.append(attr.current_time)
                else:
                    pass
            else:
                pass
        return ct

    @current_time.setter
    def current_time(self, ct):
        """Set current time for all valid attributes."""
        attr_names = dir(self)
        for attr_name in attr_names:
            if attr_name != 'current_time':
                attr = getattr(self, attr_name)

                if hasattr(attr, '_current_time_'):

                    attr.current_time = ct

                else:
                    pass

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

    def ___PRIVATE_generate_random_valid_time_instances___(self, amount=None):
        """We will generate some random valid time instances and put them in a 1d array. They can
        be used for purposes like self-checking.

        For valid time:
            - None                             : It can be everything and be changed whenever you want.
            - 'valid_only_at_its_first_instant': as it says...
            - int or float                     : Can only be this particular time instance.

        :param amount: {None, positive int}
            How many random time instances you want?  (will be functional only when it is applicable.)
        :return:
        """
        vt = self.valid_time
        if vt is None:
            if amount is None: amount = random.randint(2, 5)
            rTIs = np.random.rand(amount) * 10
        elif vt == 'valid_only_at_its_first_instant':
            rTIs = np.random.rand(1) * 10
        elif isinstance(vt, (int, float)):
            rTIs = np.array([vt,])
        else:
            raise NotImplementedError(f"valid_time = {vt} is not understandable!")

        return rTIs


    def ___Pr_compute_Ln_norm_of___(self, what, time=0, n=2, quad_degree=None):
        """We compute the :math:`L^n`-norm of the attribute name ``what`` at `time`.

        :param str what:
        :param float time:
        :param int n: We compute :math:`L^n`-norm.
        :param quad_degree:
        :return:
        """
        what = getattr(self, what)

        assert self._mesh_ is not None, " <MS> : to compute L2_norm, I need a mesh."
        if quad_degree is None:
            quad_degree = int(np.ceil((500000 / self._mesh_.elements.global_num) ** (1 / 3)))
            if quad_degree > 20: quad_degree = 20
            quad_degree = (quad_degree, quad_degree, quad_degree)

        if what.ftype == 'standard':

            _Quadrature_ = Quadrature(quad_degree, category='Gauss')
            quad_nodes, quad_weights = _Quadrature_.quad
            quad_nodes = np.meshgrid(*quad_nodes, indexing='ij')

            _, v = what.reconstruct(*quad_nodes, time=time)
            detJ = self.mesh.elements.coordinate_transformation.Jacobian(*quad_nodes)
            LOCAL = 0
            for i in v:
                local = np.sum([vij**n for vij in v[i]], axis=0)
                LOCAL += np.einsum('ijk, ijk, i, j, k -> ', local, detJ[i], *quad_weights, optimize='optimal')
            GLOBAL = COMM.allreduce(LOCAL, op=MPI.SUM) ** (1 / n)
            return GLOBAL

        else:
            raise Exception()