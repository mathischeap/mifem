# -*- coding: utf-8 -*-
"""
INTRO

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""

from components.freeze.main import FrozenOnly
from importlib import import_module
from objects.CSCG._3d.fields.vector.do.inner_product.helpers.ip1 import ___VF_INNER_PRODUCT_HELPER_1___


class _3dCSCG_Vector_Do_IP(FrozenOnly):
    def __init__(self, vf):
        self._vf_ = vf
        self._freeze_self_()

    def __call__(self, vector):
        """We compute self inner product with another vector: <self, vector>."""

        if vector.__class__.__name__ == '_3dCSCG_VectorField':

            if self._vf_.ftype == 'standard':
                if vector.ftype == 'standard':
                    w0, w1, w2 = self._vf_.func
                    u0, u1, u2 = vector.func
                    IP = ___VF_INNER_PRODUCT_HELPER_1___(w0, w1, w2, u0, u1, u2)

                    base_path = '.'.join(str(self.__class__).split(' ')[1][1:-2].split('.')[:-5]) + '.'

                    scalar_class = getattr(import_module(base_path + 'scalar.main'),
                                           '_3dCSCG_ScalarField')

                    cp_scalar = scalar_class(
                        self._vf_.mesh,
                        IP,
                        ftype='standard',
                        valid_time=self._vf_.valid_time,
                        name=self._vf_.standard_properties.name
                        + '--inner-product--'
                        + vector.standard_properties.name,
                    )
                    return cp_scalar

                else:
                    raise NotImplementedError(
                        f"a standard _3dCSCG_VectorField cannot do inner product with "
                        f"a _3dCSCG_VectorField of ftype {vector.ftype}.")
            else:
                raise NotImplementedError(f"a _3dCSCG_VectorField of ftype {self._vf_.ftype} "
                                          f"cannot do inner product.")
        else:
            raise NotImplementedError(f"a _3dCSCG_VectorField can not inner product a {vector}.")
