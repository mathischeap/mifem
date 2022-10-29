# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/10/08 8:08 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
import random
import numpy as np


class miUsTriangle_DO(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._freeze_self_()

    def generate_random_coordinates(self):
        """ We will generate some random coordinates with in this domain in a format of local mesh
        element wise. They are mainly for testing purposes.

        :return: two 1d array representing x, y coordinates.
        """
        amount = random.randint(10, 100) # 3*amount points in #amount local mesh elements.

        assert isinstance(amount, int) and amount > 0, \
            f"amount={amount} ({amount.__class__.__name__}) is wrong."

        mesh = self._mesh_
        num_local_mesh_elements = mesh.elements.num.cells

        if num_local_mesh_elements == 0:
            return np.array([]), np.array([])

        elif num_local_mesh_elements <= amount: # we will use all local mesh elements
            USED = mesh.elements.indices

        else:
            USED = random.sample(mesh.elements.indices, amount)

        xi, et = [np.random.rand(3) * 2 - 1 for _ in range(2)]
        X = list()
        Y = list()
        for i in USED:
            element = mesh.elements[i]
            x, y = element.coordinate_transformation.mapping(xi, et)
            X.append(x)
            Y.append(y)
        X = np.concatenate(X)
        Y = np.concatenate(Y)

        return X, Y














if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
