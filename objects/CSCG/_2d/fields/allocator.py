from screws.freeze.base import FrozenOnly


class _2dCSCG_FieldsAllocator(FrozenOnly):
    """"""


    @classmethod
    def ___form_name___(cls):
        return {'scalar': "_2dCSCG_ScalarField",
                'vector': "_2dCSCG_VectorField",
                }

    @classmethod
    def ___form_path___(cls):
        base_path = '.'.join(str(cls).split(' ')[1][1:-2].split('.')[:-2]) + '.'
        return {'scalar': base_path + "scalar.main",
                'vector': base_path + "vector.main",
                }