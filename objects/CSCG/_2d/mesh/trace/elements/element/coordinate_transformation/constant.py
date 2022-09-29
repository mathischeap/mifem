# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/08/29 9:20 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class _2dCSCG_TE_CT_Constant(FrozenOnly):
    """"""

    def __init__(self, te):
        """"""
        self._te_ = te
        self._Jacobian_ = True
        self._freeze_self_()

    @property
    def Jacobian(self):
        """If this trace element has a constant Jacobian, return it (a real number), otherwise
        return None.
        """
        if self._Jacobian_ is True:

            e = self._te_.CHARACTERISTIC_element
            p = self._te_.CHARACTERISTIC_position[-1]

            element = self._te_._mesh_.elements[e]

            eMark = element.type_wrt_metric.mark

            if isinstance(eMark, str) and eMark[:4] == 'Orth':
                # orthogonal(straight and parallel to axis)
                if 'x' in eMark: # not a square
                    xy = eMark[6:]
                    x, y = xy.split('y')
                    if p in 'UD':
                        LEN = float(y)
                    else:
                        LEN = float(x)
                else:
                    LEN = float(eMark[5:])

                self._Jacobian_ = LEN / 2

            # + elif, because when it is not orthogonal, the unit_normal_vector can also be constant.

            else:
                self._Jacobian_ = None


        return self._Jacobian_


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
