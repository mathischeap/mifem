from screws.freeze.base import FrozenOnly



class _3dCSCG_ADF_Allocator(FrozenOnly):
    """"""
    @classmethod
    def ___forms_name___(cls):
        return {'0-adf': "_3dCSCG_S0_ADF",
                '1-adf': "_3dCSCG_S1_ADF",
                '2-adf': "_3dCSCG_S2_ADF",
                '3-adf': "_3dCSCG_S3_ADF",

                '0-adt': "_3dCSCG_T0_ADF",
                '1-adt': "_3dCSCG_T1_ADF",
                '2-adt': "_3dCSCG_T2_ADF",

                }

    @classmethod
    def ___forms_path___(cls):
        base_path = '.'.join(str(cls).split(' ')[1][1:-2].split('.')[:-2]) + '.'

        return {'0-adf': base_path + "standard._0_AD_form",
                '1-adf': base_path + "standard._1_AD_form",
                '2-adf': base_path + "standard._2_AD_form",
                '3-adf': base_path + "standard._3_AD_form",

                '0-adt': base_path + "trace._0_AD_trace",
                '1-adt': base_path + "trace._1_AD_trace",
                '2-adt': base_path + "trace._2_AD_trace",

                }