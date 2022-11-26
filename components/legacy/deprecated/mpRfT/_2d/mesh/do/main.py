# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 
"""
import sys

if './' not in sys.path: sys.path.append('/')

import random
from components.freeze.base import FrozenOnly
from objects.mpRfT._2d.mesh.do.find import mpRfT2_Mesh_Do_Find
import numpy as np

class mpRfT2_Mesh_Do(FrozenOnly):
    """"""
    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._find_ = None
        self._freeze_self_()

    @property
    def find(self):
        """"""
        if self._find_ is None:
            self._find_ = mpRfT2_Mesh_Do_Find(self._mesh_)
        return self._find_

    def evolve(self):
        """Applying all the refinements to the cscg mesh to make a new mpRfT2 mesh."""
        rfd = self._mesh_.refinements.future.rfd
        return self._mesh_.__class__(self._mesh_.cscg, self._mesh_.dN, rfd)




    def generate_random_coordinates(self):
        """ We will generate some random coordinates with in this domain in a format of local mesh
        element wise. They are mainly for testing purposes.

        :return: two 1d array representing x, y coordinates.
        """
        amount = random.randint(10, 100) # 3*amount points in #amount local mesh elements.

        assert isinstance(amount, int) and amount > 0, \
            f"amount={amount} ({amount.__class__.__name__}) is wrong."

        mesh = self._mesh_
        num_local_rc = len(mesh)

        if num_local_rc == 0:
            return np.array([]), np.array([])

        elif num_local_rc <= amount: # we will use all local mesh elements
            USED = mesh.rcfc

        else:
            USED = random.sample(mesh.rcfc.keys(), amount)

        xi, et = [np.random.rand(3) * 2 - 1 for _ in range(2)]
        X = list()
        Y = list()
        for rc_rp in USED:
            root_cell = mesh[rc_rp]
            x, y = root_cell.coordinate_transformation.mapping(xi, et)
            X.append(x)
            Y.append(y)
        X = np.concatenate(X)
        Y = np.concatenate(Y)

        return X, Y

if __name__ == '__main__':
    # mpiexec -n 4 python objects/mpRfT/_2d/mesh/do/main.py
    pass