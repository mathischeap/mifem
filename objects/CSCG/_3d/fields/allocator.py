# -*- coding: utf-8 -*-
from screws.freeze.base import FrozenOnly



class _3dCSCG_Field_Allocator(FrozenOnly):
    """"""
    @classmethod
    def ___forms_name___(cls):
        return {'scalar': "_3dCSCG_ScalarField",
                'vector': "_3dCSCG_VectorField",
                'tensor': "_3dCSCG_TensorField",
                }

    @classmethod
    def ___forms_path___(cls):
        base_path = '.'.join(str(cls).split(' ')[1][1:-2].split('.')[:-2]) + '.'

        return {'scalar': base_path + "scalar.main",
                'vector': base_path + "vector.main",
                'tensor': base_path + "tensor.main",
                }