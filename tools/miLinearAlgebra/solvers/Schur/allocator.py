# -*- coding: utf-8 -*-
"""
A collection of all parallel solvers.
"""
from importlib import import_module

from components.freeze.main import FrozenOnly


class SchurSolverDistributor(FrozenOnly):
    """"""
    def __init__(self, rank, blocks, routine='auto', name='SchurSolver'):
        """"""
        solver_name = self.___solver_name___()[rank]
        solver_path = self.___solver_path___()[rank]
        self._solver_ = getattr(import_module(solver_path), solver_name)(rank, blocks, routine, name)

    def __call__(self, A, b,  *args, **kwargs):
        """
        :return: Return a tuple of 5 outputs:

            1. (LocallyFullVector) results -- The result vector.
            2. (int) info -- The info which provides convergence information:

                * 0 : successful exit
                * >0 : convergence to tolerance not achieved, number of iterations
                * -1 : divergence

            3. (float) beta -- The residual.
            4. (int) ITER -- The number of outer iterations.
            5. (str) message

        """
        assert A.__class__.__name__ == "EWC_SparseMatrix", \
            f"A needs to be a EWC_SparseMatrix'. Now I get {A.__class__}."
        assert b.__class__.__name__ == "EWC_ColumnVector", \
            f"b needs to be a 'EWC_ColumnVector'. Now I get {b.__class__}."
        assert A.bmat_shape is not False and b.con_shape is not False
        Results = self._solver_(A, b, *args, **kwargs)
        assert Results[0].__class__.__name__ == 'dict', \
            f"results must be in a dict!"
        return Results

    @classmethod
    def ___solver_name___(cls):
        """"""
        return {2: 'Rank2', }

    @classmethod
    def ___solver_path___(cls):
        base_path = '.'.join(str(cls).split(' ')[1][1:-2].split('.')[:-2]) + '.'
        return {
            2: base_path + 'rank2.main',
        }
