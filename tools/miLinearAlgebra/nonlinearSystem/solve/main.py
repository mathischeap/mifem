# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/7/24 16:28
"""
from components.freeze.main import FrozenOnly

from tools.miLinearAlgebra.nonlinearSystem.solve.NewtonRaphson.main import nLS_Solve_NewtonRaphson


class NonLinearSystem_Solve(FrozenOnly):
    def __init__(self, nLS):
        """"""
        self._nLS_ = nLS
        self._nLS_solver_ = 'Newton-Raphson'
        self._NewtonRaphson_ = nLS_Solve_NewtonRaphson(nLS)
        self._freeze_self_()

    @property
    def solver(self):
        return self._nLS_solver_

    @solver.setter
    def solver(self, solver):
        assert isinstance(solver, str), f"I need a str solver name."
        self._nLS_solver_ = solver

    @property
    def Newton_Raphson(self):
        """The Newton-Raphson solver instance."""
        return self._NewtonRaphson_

    def __call__(self, *args, **kwargs):
        """"""
        if self._nLS_solver_ == 'Newton-Raphson':
            RES = self._NewtonRaphson_(*args, **kwargs)
        else:
            raise NotImplementedError()

        return RES
