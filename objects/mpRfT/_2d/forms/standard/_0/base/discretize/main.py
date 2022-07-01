# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/16 4:14 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.mpRfT._2d.forms.standard._0.base.discretize.scalar.standard import mpRfT2_S0F_Discretize_Standard_Scalar



class mpRfT2_S0F_Discretize(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._Pr_standard_scalar = mpRfT2_S0F_Discretize_Standard_Scalar(f)
        self._freeze_self_()

    def __call__(self, target='analytic_expression'):
        """

        Parameters
        ----------
        target : str

        Returns
        -------

        """

        if target == 'analytic_expression':

            if self._f_.analytic_expression.__class__.__name__ == 'mpRfT2_Scalar':

                if self._f_.analytic_expression.ftype == 'standard':
                    LCC =  self._Pr_standard_scalar(target)
                else:
                    raise NotImplementedError(f"mpRfT2_S0F cannot (target func) "
                                              f"discretize mpRfT2_S0F of ftype={self._f_.analytic_expression.ftype}")

            else:
                raise NotImplementedError(f'mpRfT2_S0F can not (target func) '
                                          f'discretize {self._f_.analytic_expression.__class__}.')

        else:
            raise NotImplementedError(f"mpRfT2_S0F cannot discretize "
                                      f"while targeting at {target}.")

        self._f_.cochain.local = LCC




if __name__ == "__main__":
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/form/standard/_0/base/discretize/main.py
    pass
