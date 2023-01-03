# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/8/4 21:15
"""
from components.freeze.main import FrozenOnly
from tools.miLinearAlgebra.dataStructures.vectors.locallyFull.main import LocallyFullVector
from tools.miLinearAlgebra.nonlinearSystem.solve.NewtonRaphson.regular import nLS_Solve_NR_regular


class nLS_Solve_NewtonRaphson(FrozenOnly):
    def __init__(self, nLS):
        """"""
        self._nLS_ = nLS
        self._routine_ = 'regular'
        self._regular_ = nLS_Solve_NR_regular(nLS)
        self._freeze_self_()

    @property
    def routine(self):
        """The routine this Newton-Raphson scheme is using."""
        return self._routine_

    @routine.setter
    def routine(self, routine):
        """"""
        self._routine_ = routine

    def __call__(self, x0, **kwargs):
        """ """
        x0 = self.___PRIVATE_parse_x0___(x0)

        if self.routine == 'regular':
            # we call the regular solve to solve it.
            RES = self._regular_(x0, **kwargs)
        else:
            raise NotImplementedError(f"routine={self.routine} is not implemented.")

        return RES

    def ___PRIVATE_parse_x0___(self, x0):

        if x0 == 0:  # we make it an empty LocallyFullVector
            x0 = LocallyFullVector(self._nLS_.num.equations)
        else:
            assert x0.__class__.__name__ == 'LocallyFullVector', \
                f"x0={x0} is wrong, I need a LocallyFullVector."

        return x0
