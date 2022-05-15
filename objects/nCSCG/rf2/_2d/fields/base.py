# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/13 2:55 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.nCSCG.rf2.base.fields.base import nCSCG_FieldBase


class _2nCSCG_FieldBase(nCSCG_FieldBase):
    """"""
    def __init__(self, mesh, ftype, valid_time, name):
        """"""
        super(_2nCSCG_FieldBase, self).__init__(mesh, ftype=ftype, valid_time=valid_time, name=name)
        self.standard_properties.___PRIVATE_add_tag___('_2nCSCG_RF2_Field')





if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
