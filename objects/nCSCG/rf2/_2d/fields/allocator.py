



from screws.freeze.base import FrozenOnly


class _2nCSCG_RF2_FieldsAllocator(FrozenOnly):
    """"""


    @classmethod
    def ___form_name___(cls):
        return {'scalar': "_2nCSCG_RF2_ScalarField",
                'vector': "_2nCSCG_RF2_VectorField",
                }

    @classmethod
    def ___form_path___(cls):
        base_path = '.'.join(str(cls).split(' ')[1][1:-2].split('.')[:-2]) + '.'
        return {'scalar': base_path + "scalar.main",
                'vector': base_path + "vector.main",
                }