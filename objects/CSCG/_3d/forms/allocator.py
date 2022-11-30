# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly



class _3dCSCG_SF_Allocator(FrozenOnly):
    """"""

    @classmethod
    def ___forms_name___(cls):
        return {
            '3-f': "_3dCSCG_3Form",
            '2-f': "_3dCSCG_2Form",
            '1-f': "_3dCSCG_1Form",
            '0-f': "_3dCSCG_0Form",

            '0-t': "_3dCSCG_0Trace",
            '1-t': "_3dCSCG_1Trace",
            '2-t': "_3dCSCG_2Trace",

            '0-e': "_3dCSCG_0Edge",
            '1-e': "_3dCSCG_1Edge",

            '0-lt': "_3dCSCG_0LocalTrace",
            '1-lt': "_3dCSCG_1LocalTrace",
            '2-lt': "_3dCSCG_2LocalTrace",
        }

    @classmethod
    def ___forms_path___(cls):
        base_path = '.'.join(str(cls).split(' ')[1][1:-2].split('.')[:-2]) + '.'

        return {
            '3-f': base_path + "standard._3s.main",
            '2-f': base_path + "standard._2s.main",
            '1-f': base_path + "standard._1s.main",
            '0-f': base_path + "standard._0s.main",

            '0-t': base_path + "trace._0tr.main",
            '1-t': base_path + "trace._1tr.main",
            '2-t': base_path + "trace._2tr.main",

            '0-e': base_path + "edge._0eg.main",
            '1-e': base_path + "edge._1eg.main",

            '0-lt': base_path + "localTrace._0ltf.main",
            '1-lt': base_path + "localTrace._1ltf.main",
            '2-lt': base_path + "localTrace._2ltf.main",
        }