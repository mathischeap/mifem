# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/19 3:31 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly


class mpRfT2_FormsAllocator(FrozenOnly):
    """"""


    @classmethod
    def ___form_name___(cls):
        return {'0-f-i': "mpRfT2_Si0F", # d on inner 0-form is grad
                '0-f-o': "mpRfT2_So0F", # d on outer 0-form is curl

                '1-f-i': "mpRfT2_Si1F",  # d on inner 1-form is rot
                '1-f-o': "mpRfT2_So1F",  # d on outer 1-form is div

                '2-f-i': "mpRfT2_Si2F",
                '2-f-o': "mpRfT2_So2F",

                'nst': "mpRfT2_NSgF",
                'est': "mpRfT2_ESgF",
                }

    @classmethod
    def ___form_path___(cls):
        base_path = '.'.join(str(cls).split(' ')[1][1:-2].split('.')[:-2]) + '.'
        return {'0-f-i': base_path + "standard._0.inner.main", # d on inner 0-form is grad
                '0-f-o': base_path + "standard._0.outer.main", # d on outer 0-form is curl

                '1-f-i': base_path + "standard._1.inner.main",  # d on inner 1-form is rot
                '1-f-o': base_path + "standard._1.outer.main",  # d on outer 1-form is div

                '2-f-i': base_path + "standard._2.inner.main",
                '2-f-o': base_path + "standard._2.outer.main",

                'nst': base_path + "segment.node.main",
                'est': base_path + "segment.edge.main",
                }


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
