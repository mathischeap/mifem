# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 10/5/2022 2:49 PM
"""
from components.freeze.main import FrozenOnly

from objects.CSCG._2d.fields.scalar.do.cross_product.helpers.helper2 import ___SF_CROSS_PRODUCT_HELPER_2___
from objects.CSCG._2d.fields.scalar.do.cross_product.helpers.helper1 import ___SF_CROSS_PRODUCT_HELPER_1___

from importlib import import_module


class miUsTriangle_ScalarField_CrossProduct(FrozenOnly):
    """w = [0 0 w]^T, U = [u, v, 0]^T, A = w X U = [-wv wu 0]^T"""
    def __init__(self, sf):
        self._sf_ = sf
        self._freeze_self_()

    def __call__(self, vector):
        """We compute self cross_product another vector: self(a scalar) x vector"""

        # w = [0 0 w]^T, U = [u, v, 0]^T, A = w X U = [-wv wu 0]^T
        if vector.__class__.__name__ == 'miUsGrid_Triangular_Vector':

            if self._sf_.ftype == 'standard':
                if vector.ftype == 'standard':
                    w = self._sf_.func[0]
                    u, v = vector.func
                    CP0 = ___SF_CROSS_PRODUCT_HELPER_1___(w, v)
                    CP1 = ___SF_CROSS_PRODUCT_HELPER_2___(w, u)

                    base_path = '.'.join(str(self.__class__).split(' ')[1][1:-2].split('.')[:-4]) + '.'

                    vector_class = getattr(import_module(base_path + 'vector.main'),
                                           'miUsGrid_Triangular_Vector')

                    cp_vector = vector_class(
                        self._sf_.mesh,
                        [CP0, CP1],
                        ftype='standard',
                        valid_time=self._sf_.valid_time,
                        name=self._sf_.standard_properties.name +
                              '--cross-X-product--' +
                              vector.standard_properties.name
                    )
                    return cp_vector
                else:
                    raise NotImplementedError(
                        f"a standard miUsGrid_Triangular_Scalar cannot do cross product with "
                        f"a miUsGrid_Triangular_Vector of ftype {vector.ftype}.")
            else:
                raise NotImplementedError(f"a miUsGrid_Triangular_Scalar of ftype {self._sf_.ftype} "
                                          f"cannot do cross product.")
        else:
            raise NotImplementedError(f"a miUsGrid_Triangular_Vector can not cross product a {vector}.")
