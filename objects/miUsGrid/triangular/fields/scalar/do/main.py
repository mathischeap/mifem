# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 2:55 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.miUsGrid.triangular.fields.scalar.do.cross_product import miUsTriangle_ScalarField_CrossProduct
from objects.miUsGrid.triangular.fields.scalar.do.inner_product import miUsTriangle_ScalarField_InnerProduct


class miUsGrid_Triangular_Scalar_Do(FrozenOnly):
    """A collection of other methods beside the standard ones like reconstruct, discretize, visualize, etc."""

    def __init__(self, scalar):
        """"""
        self._scalar_ = scalar
        self._cp_ = miUsTriangle_ScalarField_CrossProduct(scalar)
        self._ip_ = miUsTriangle_ScalarField_InnerProduct(scalar)
        self._freeze_self_()


    def evaluate_func_at_time(self, time=None):
        return self._scalar_.___DO_evaluate_func_at_time___(time=time)

    @property
    def cross_product(self):
        return self._cp_

    @property
    def inner_product(self):
        return self._ip_


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
