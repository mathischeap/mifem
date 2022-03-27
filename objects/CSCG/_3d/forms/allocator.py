from screws.freeze.base import FrozenOnly



class _3dCSCG_SF_Allocator(FrozenOnly):
    """"""
    @classmethod
    def ___forms_name___(cls):
        return {'3-f': "_3dCSCG_3Form",
                '2-f': "_3dCSCG_2Form",
                '1-f': "_3dCSCG_1Form",
                '0-f': "_3dCSCG_0Form",

                '0-t': "_3dCSCG_0Trace",
                '1-t': "_3dCSCG_1Trace",
                '2-t': "_3dCSCG_2Trace",

                '0-e': "_3dCSCG_0Edge",
                '1-e': "_3dCSCG_1Edge",

                '0-n': "_3dCSCG_0Node",
                }

    @classmethod
    def ___forms_path___(cls):
        base_path = '.'.join(str(cls).split(' ')[1][1:-2].split('.')[:-2]) + '.'

        return {'3-f': base_path + "standard._3_form.main",
                '2-f': base_path + "standard._2_form.main",
                '1-f': base_path + "standard._1_form.main",
                '0-f': base_path + "standard._0_form.main",

                '0-t': base_path + "trace._0_trace.main",
                '1-t': base_path + "trace._1_trace.main",
                '2-t': base_path + "trace._2_trace.main",

                '0-e': base_path + "edge._0_edge",
                '1-e': base_path + "edge._1_edge",

                '0-n': base_path + "node._0_node",
                }