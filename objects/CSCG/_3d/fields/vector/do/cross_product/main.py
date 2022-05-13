# -*- coding: utf-8 -*-
from screws.freeze.main import FrozenOnly
from importlib import import_module
from objects.CSCG._3d.fields.vector.do.cross_product.helpers.cp1 import ___VF_CROSS_PRODUCT_HELPER_1___




class _3dCSCG_Vector_Do_CP(FrozenOnly):
    def __init__(self, vf):
        self._vf_ = vf
        self._freeze_self_()

    def __call__(self, vector):
        """We compute self cross_product another vector: self X vector"""

        if vector.__class__.__name__ == '_3dCSCG_VectorField':

            if self._vf_.ftype == 'standard':
                if vector.ftype == 'standard':
                    w0, w1, w2 = self._vf_.func
                    u0, u1, u2 = vector.func
                    CP0 = ___VF_CROSS_PRODUCT_HELPER_1___(w1, u2, w2, u1)
                    CP1 = ___VF_CROSS_PRODUCT_HELPER_1___(w2, u0, w0, u2)
                    CP2 = ___VF_CROSS_PRODUCT_HELPER_1___(w0, u1, w1, u0)

                    base_path = '.'.join(str(self.__class__).split(' ')[1][1:-2].split('.')[:-5]) + '.'

                    vector_class = getattr(import_module(base_path + 'vector.main'),
                                           '_3dCSCG_VectorField')

                    cp_vector = vector_class(self._vf_.mesh,
                                                    [CP0, CP1, CP2],
                                                    ftype='standard',
                                                    valid_time=self._vf_.valid_time,
                                                    name = self._vf_.standard_properties.name
                                                         + '--cross-X-product--'
                                                         + vector.standard_properties.name
                                                    )
                    return cp_vector
                else:
                    raise NotImplementedError(
                        f"a standard _3dCSCG_VectorField cannot do cross product with "
                        f"a _3dCSCG_VectorField of ftype {vector.ftype}.")
            else:
                raise NotImplementedError(f"a _3dCSCG_VectorField of ftype {self._vf_.ftype} "
                                          f"cannot do cross product.")
        else:
            raise NotImplementedError(f"a _3dCSCG_VectorField can not cross product a {vector}.")

