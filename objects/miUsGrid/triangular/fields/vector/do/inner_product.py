# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 10/5/2022 2:49 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from screws.freeze.main import FrozenOnly

from objects.CSCG._2d.fields.scalar.do.cross_product.helpers.helper2 import ___SF_CROSS_PRODUCT_HELPER_2___

from importlib import import_module

class miUsTriangle_VectorField_InnerProduct(FrozenOnly):
    def __init__(self, vf):
        self._vf_ = vf
        self._freeze_self_()

    def __call__(self, vector):
        if vector.__class__.__name__ == 'miUsGrid_Triangular_Vector':

            if self._vf_.ftype == 'standard':
                if vector.ftype == 'standard':
                    u, v = self._vf_.func
                    a, b = vector.func
                    IP1 = ___SF_CROSS_PRODUCT_HELPER_2___(u, a)
                    IP2 = ___SF_CROSS_PRODUCT_HELPER_2___(v, b)

                    base_path = '.'.join(str(self.__class__).split(' ')[1][1:-2].split('.')[:-4]) + '.'

                    vector_class = getattr(import_module(base_path + 'vector.main'),
                                           'miUsGrid_Triangular_Vector')

                    Ip_vector = vector_class(
                        self._vf_.mesh,
                        [IP1, IP2],
                        ftype='standard',
                        valid_time=self._vf_.valid_time,
                        name =self._vf_.standard_properties.name
                              + '--inner-product--'
                              + vector.standard_properties.name
                    )
                    return Ip_vector
                else:
                    raise NotImplementedError(
                        f"a standard miUsGrid_Triangular_Vector cannot do inner product with "
                        f"a miUsGrid_Triangular_Vector of ftype {vector.ftype}.")
            else:
                raise NotImplementedError(f"a miUsGrid_Triangular_Vector of ftype {self._vf_.ftype} "
                                          f"cannot do inner product.")
        else:
            raise NotImplementedError(f"a miUsGrid_Triangular_Vector can not inner product a {vector}.")


if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
