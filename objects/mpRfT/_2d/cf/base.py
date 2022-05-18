# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/13 2:55 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.mpRfT.base.cf.base import mpRfT_ContinuousField_Base


class mpRfT2_ContinuousField(mpRfT_ContinuousField_Base):
    """"""
    def __init__(self, mesh, ftype, valid_time, name):
        """"""
        assert mesh.__class__.__name__ == 'mpRfT2_Mesh', f"mesh={mesh} is wrong."
        super(mpRfT2_ContinuousField, self).__init__(mesh, ftype=ftype, valid_time=valid_time, name=name)






if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
