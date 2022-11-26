# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly



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

        return {'0-adf': base_path + "standard._0s.main",
                '1-adf': base_path + "standard._1s.main",
                '2-adf': base_path + "standard._2s.main",
                '3-adf': base_path + "standard._3s.main",

                '0-adt': base_path + "trace._0tr.main",
                '1-adt': base_path + "trace._1tr.main",
                '2-adt': base_path + "trace._2tr.main",
                }