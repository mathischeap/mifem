# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/20 2:36 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from components.freeze.base import FrozenOnly
from root.config.main import RANK, MASTER_RANK

class miUsGrid_FormsAllocator(FrozenOnly):
    """"""

    @classmethod
    def listing(cls, printing=True, returning=True):
        """"""
        if RANK != MASTER_RANK: return

        listing = ''
        names = cls.___form_name___()
        paths = cls.___form_path___()
        for ID in names:
            assert ID in paths, f"ID={ID} finds no path."
            listing += ">>> " + ID + ' ~ ' + names[ID] + '\n\n'

        if printing:
            print(listing)
        else:
            pass

        if returning:
            return listing
        else:
            pass

    @classmethod
    def ___form_name___(cls):
        return {'0-f-i': "miUsTriangular_S0F_Inner", # d on inner 0-form is grad
                '1-f-i': "miUsTriangular_S1F_Inner", # d on inner 1-form is rot
                '2-f-i': "miUsTriangular_S2F_Inner",

                '0-f-o': "miUsTriangular_S0F_Outer", # d on outer 0-form is curl
                '1-f-o': "miUsTriangular_S1F_Outer", # d on outer 1-form is div
                '2-f-o': "miUsTriangular_S2F_Outer",
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



if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
