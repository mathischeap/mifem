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

class miUsTriangle_ScalarField_InnerProduct(FrozenOnly):
    """w = [0 0 w]^T, U = [u, v, 0]^T, A = w X U = [-wv wu 0]^T"""
    def __init__(self, sf):
        self._sf_ = sf
        self._freeze_self_()

    def __call__(self, scalar):
        """We compute self inner product another scalar: self(a scalar) cdot vector"""

        # w = w, U = u, A = w cdot U = wu
        if scalar.__class__.__name__ == 'miUsGrid_Triangular_Scalar':

            if self._sf_.ftype == 'standard':
                if scalar.ftype == 'standard':
                    w = self._sf_.func[0]
                    u = scalar.func[0]
                    IP = ___SF_CROSS_PRODUCT_HELPER_2___(w, u)

                    base_path = '.'.join(str(self.__class__).split(' ')[1][1:-2].split('.')[:-4]) + '.'

                    vector_class = getattr(import_module(base_path + 'scalar.main'),
                                           'miUsGrid_Triangular_Scalar')

                    cp_vector = vector_class(self._sf_.mesh,
                                                    IP,
                                                    ftype='standard',
                                                    valid_time=self._sf_.valid_time,
                                                    name =self._sf_.standard_properties.name
                                                          + '--inner-product--'
                                                          + scalar.standard_properties.name
                                                    )
                    return cp_vector
                else:
                    raise NotImplementedError(
                        f"a standard miUsGrid_Triangular_Scalar cannot do inner product with "
                        f"a miUsGrid_Triangular_Scalar of ftype {scalar.ftype}.")

            else:
                raise NotImplementedError(f"a miUsGrid_Triangular_Scalar of ftype {self._sf_.ftype} "
                                          f"cannot do inner product.")
        else:
            raise NotImplementedError(f"a miUsGrid_Triangular_Scalar can not inner product a {scalar}.")


if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
