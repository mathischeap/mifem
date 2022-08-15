# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/8/11 0:36
"""
import sys

if './' not in sys.path: sys.path.append('./')
from screws.freeze.main import FrozenOnly
from importlib import import_module

from objects.CSCG._2d.fields.vector.do.inner_product.helpers.helper import VF2_IP_H


class _2CSCG_VectorField_InnerProduct(FrozenOnly):
    def __init__(self, vf):
        """"""
        self._vf_ = vf
        self._freeze_self_()

    def __call__(self, vector):
        """We compute self cross_product another vector: self(a scalar) x vector"""

        if vector.__class__.__name__ == '_2dCSCG_VectorField':

            if self._vf_.ftype == 'standard':
                if vector.ftype == 'standard':
                    f0, f1 = self._vf_.func
                    f2, f3 = vector.func
                    IP = VF2_IP_H(f0, f1, f2, f3)

                    base_path = '.'.join(str(self.__class__).split(' ')[1][1:-2].split('.')[:-5]) + '.'

                    vector_class = getattr(import_module(base_path + 'scalar.main'),
                                           '_2dCSCG_ScalarField')

                    cp_vector = vector_class(self._vf_.mesh,
                                                    IP,
                                                    ftype='standard',
                                                    valid_time=self._vf_.valid_time,
                                                    name = self._vf_.standard_properties.name
                                                         + '--inner-product--'
                                                         + vector.standard_properties.name
                                                    )
                    return cp_vector
                else:
                    raise NotImplementedError(
                        f"a standard _2dCSCG_VectorField cannot do inner product with "
                        f"a _2dCSCG_VectorField of ftype {vector.ftype}.")
            else:
                raise NotImplementedError(f"a _2dCSCG_VectorField of ftype {self._vf_.ftype} "
                                          f"cannot do inner product.")
        else:
            raise NotImplementedError(f"a _2dCSCG_VectorField can not inner product a {vector}.")



if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
