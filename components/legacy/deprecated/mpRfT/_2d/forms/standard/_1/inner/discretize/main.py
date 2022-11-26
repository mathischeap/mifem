# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/24 12:26 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly
from objects.mpRfT._2d.forms.standard._1.inner.discretize.vector.standard import mpRfT2_Si1F_Discretize_Standard_Vector


class mpRfT2_Si1F_Discretize(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._Pr_standard_vector = mpRfT2_Si1F_Discretize_Standard_Vector(f)
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

            if self._f_.analytic_expression.__class__.__name__ == 'mpRfT2_Vector':

                if self._f_.analytic_expression.ftype == 'standard':
                    LCC =  self._Pr_standard_vector(target)
                else:
                    raise NotImplementedError(f"mpRfT2_Si1F cannot (target analytic_expression) "
                                              f"discretize mpRfT2_Si1F of ftype={self._f_.analytic_expression.ftype}")

            else:
                raise NotImplementedError(f'mpRfT2_Si1F can not (target analytic_expression) '
                                          f'discretize {self._f_.analytic_expression.__class__}.')

        else:
            raise NotImplementedError(f"mpRfT2_Si1F cannot discretize "
                                      f"while targeting at {target}.")

        self._f_.cochain.local = LCC



if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
