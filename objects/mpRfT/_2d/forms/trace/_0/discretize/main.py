# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/16 4:14 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.mpRfT._2d.forms.trace._0.discretize.scalar.standard import mpRfT2_T0F_Discretize_Standard_Scalar



class mpRfT2_T0F_Discretize(FrozenOnly):
    """"""

    def __init__(self, t):
        """"""
        self._t_ = t
        self._Pr_standard_scalar = mpRfT2_T0F_Discretize_Standard_Scalar(t)
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

            if self._t_.TW.func.__class__.__name__ == 'mpRfT2_Scalar':

                if self._t_.TW.func.ftype == 'standard':
                    LCC =  self._Pr_standard_scalar(target)
                else:
                    raise NotImplementedError(f"mpRfT2_T0F cannot (target func) "
                                              f"discretize mpRfT2_T1F of ftype={self._t_.TW.func.ftype}")

            else:
                raise NotImplementedError(f'mpRfT2_T0F can not (target func) '
                                          f'discretize {self._t_.TW.func.__class__}.')

        else:
            raise NotImplementedError(f"mpRfT2_T0F cannot discretize "
                                      f"while targeting at {target}.")

        if update_cochain: self._t_.cochain.local = LCC

        return LCC




if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/trace/_0/discretize/main.py
    from __init__ import rfT2

    fc = rfT2.rf(100)

    t0 = fc('0-t')
