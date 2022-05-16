# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class _2nCSCG_RF2_FormsAllocator(FrozenOnly):
    """"""


    @classmethod
    def ___form_name___(cls):
        return {'0-f-i': "_2nCSCG_RF2_InnerStandard0Form", # d on inner 0-form is grad
                '1-f-i': "_2nCSCG_RF2_InnerStandard1Form", # d on inner 1-form is rot
                '2-f-i': "_2nCSCG_RF2_InnerStandard2Form",

                '0-f-o': "_2nCSCG_RF2_OuterStandard0Form", # d on outer 0-form is curl
                '1-f-o': "_2nCSCG_RF2_OuterStandard1Form", # d on outer 1-form is div
                '2-f-o': "_2nCSCG_RF2_OuterStandard2Form",

                }

    @classmethod
    def ___form_path___(cls):
        base_path = '.'.join(str(cls).split(' ')[1][1:-2].split('.')[:-2]) + '.'
        return {'0-f-i': base_path + "standard._0.inner.main", # d on inner 0-form is grad
                '1-f-i': base_path + "standard._1.inner.main", # d on inner 1-form is rot
                '2-f-i': base_path + "standard._2.inner.main",

                '0-f-o': base_path + "standard._0.outer.main", # d on outer 0-form is curl
                '1-f-o': base_path + "standard._1.outer.main", # d on outer 1-form is div
                '2-f-o': base_path + "standard._2.outer.main",

                }







if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
