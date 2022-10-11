# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 2:55 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.miUsGrid.triangular.fields.vector.do.inner_product import miUsTriangle_VectorField_InnerProduct


class miUsGrid_Triangular_Vector_Do(FrozenOnly):
    """A collection of other methods beside the standard ones like reconstruct, discretize, visualize, etc."""

    def __init__(self, vector):
        """"""
        self._vector_ = vector
        self._ip_ = miUsTriangle_VectorField_InnerProduct(vector)
        self._freeze_self_()


    def evaluate_func_at_time(self, time=None):
        return self._vector_.___DO_evaluate_func_at_time___(time=time)

    @property
    def inner_product(self):
        return self._ip_


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
