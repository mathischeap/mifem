# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/16 4:14 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly
from objects.mpRfT._2d.forms.segment.node.discretize.scalar.standard import mpRfT2_NSgF_Discretize_Standard_Scalar





class mpRfT2_NSgF_Discretize(FrozenOnly):
    """"""

    def __init__(self, t):
        """"""
        self._t_ = t
        self._Pr_standard_scalar = mpRfT2_NSgF_Discretize_Standard_Scalar(t)
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

            if self._t_.analytic_expression.__class__.__name__ == 'mpRfT2_Scalar':

                if self._t_.analytic_expression.ftype == 'standard':
                    sgw_LCC =  self._Pr_standard_scalar(target)

                else:
                    raise NotImplementedError(f"mpRfT2_NSgF cannot (target analytic_expression) "
                                              f"discretize mpRfT2_Scalar of ftype="
                                              f"{self._t_.analytic_expression.ftype}")

            else:
                raise NotImplementedError(f'mpRfT2_NSgF can not (target analytic_expression) '
                                          f'discretize {self._t_.analytic_expression.__class__}.')

        elif target == 'boundary_condition':

            BC = self._t_.BC

            if BC.analytic_expression.__class__.__name__ == 'mpRfT2_Scalar':

                if BC.analytic_expression.ftype == 'standard':
                    sgw_LCC =  self._Pr_standard_scalar(target)

                else:
                    raise NotImplementedError(f"mpRfT2_NSgF cannot (target BC) "
                                              f"discretize mpRfT2_Scalar of ftype="
                                              f"{BC.analytic_expression.ftype}")

            else:
                raise NotImplementedError(f'mpRfT2_NSgF can not (target BC) '
                                          f'discretize {BC.analytic_expression.__class__}.')


        else:
            raise NotImplementedError(f"mpRfT2_NSgF cannot discretize "
                                      f"while targeting at {target}.")

        self._t_.cochain.trace = sgw_LCC







if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/trace/_0/discretize/main.py
    from __init__ import rfT2

    fc = rfT2.rf(100)

    t0 = fc('0-t')
