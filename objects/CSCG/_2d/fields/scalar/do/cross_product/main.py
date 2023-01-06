# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/8/11 0:25
"""
import sys

if './' not in sys.path:
    sys.path.append('./')
from components.freeze.main import FrozenOnly

from objects.CSCG._2d.fields.scalar.do.cross_product.helpers.helper2 import ___SF_CROSS_PRODUCT_HELPER_2___
from objects.CSCG._2d.fields.scalar.do.cross_product.helpers.helper1 import ___SF_CROSS_PRODUCT_HELPER_1___

from importlib import import_module


class _2dCSCG_SclarField_CrossProduct(FrozenOnly):
    """w = [0 0 w]^T, U = [u, v, 0]^T, A = w X U = [-wv wu 0]^T"""
    def __init__(self, sf):
        self._sf_ = sf
        self._freeze_self_()

    def __call__(self, vector):
        """We compute self cross_product another vector: self(a scalar) x vector"""

        # w = [0 0 w]^T, U = [u, v, 0]^T, A = w X U = [-wv wu 0]^T
        if vector.__class__.__name__ == '_2dCSCG_VectorField':

            if self._sf_.ftype == 'standard':
                if vector.ftype == 'standard':
                    w = self._sf_.func[0]
                    u, v = vector.func
                    CP0 = ___SF_CROSS_PRODUCT_HELPER_1___(w, v)
                    CP1 = ___SF_CROSS_PRODUCT_HELPER_2___(w, u)

                    base_path = '.'.join(str(self.__class__).split(' ')[1][1:-2].split('.')[:-5]) + '.'

                    vector_class = getattr(import_module(base_path + 'vector.main'),
                                           '_2dCSCG_VectorField')

                    cp_vector = vector_class(
                        self._sf_.mesh,
                        [CP0, CP1],
                        ftype='standard',
                        valid_time=self._sf_.valid_time,
                        name=self._sf_.standard_properties.name
                              + '--cross-X-product--'
                              + vector.standard_properties.name
                    )
                    return cp_vector
                else:
                    raise NotImplementedError(
                        f"a standard _2dCSCG_VectorField cannot do cross product with "
                        f"a _2dCSCG_VectorField of ftype {vector.ftype}.")
            else:
                raise NotImplementedError(f"a _2dCSCG_VectorField of ftype {self._sf_.ftype} "
                                          f"cannot do cross product.")
        else:
            raise NotImplementedError(f"a _2dCSCG_VectorField can not cross product a {vector}.")


if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
