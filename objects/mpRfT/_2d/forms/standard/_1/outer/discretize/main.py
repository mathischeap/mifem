# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/24 12:26 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.mpRfT._2d.forms.standard._1.outer.discretize.vector.standard import mpRfT2_So1F_Discretize_Standard_Vector


class mpRfT2_So1F_Discretize(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._Pr_standard_vector = mpRfT2_So1F_Discretize_Standard_Vector(f)
        self._freeze_self_()


    def __call__(self, update_cochain=True, target='func'):
        """

        Parameters
        ----------
        update_cochain : bool
            If we send the output to `cochain.local`.
        target : str

        Returns
        -------

        """

        if target == 'func':

            if self._f_.TW.func.__class__.__name__ == 'mpRfT2_Vector':

                if self._f_.TW.func.ftype == 'standard':
                    LCC =  self._Pr_standard_vector(target)
                else:
                    raise NotImplementedError(f"mpRfT2_So1F cannot (target func) "
                                              f"discretize mpRfT2_So1F of ftype={self._f_.TW.func.ftype}")

            else:
                raise NotImplementedError(f'mpRfT2_So1F can not (target func) '
                                          f'discretize {self._f_.TW.func.__class__}.')

        else:
            raise NotImplementedError(f"mpRfT2_So1F cannot discretize "
                                      f"while targeting at {target}.")

        if update_cochain: self._f_.cochain.local = LCC

        return LCC


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
