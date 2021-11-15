"""
A collection of all parallel solvers.
"""



import sys
if './' not in sys.path: sys.path.append('./')


from importlib import import_module


from SCREWS.frozen import FrozenOnly


class ParallelSolverDistributor(FrozenOnly):
    """"""
    def __init__(self, solver_name):
        """"""
        self._sn_ = solver_name
        solver_path = self.___parallel_solver_path___() + solver_name
        self._solver_ = getattr(import_module(solver_path), 'solve')

    def __call__(self, *args, **kwargs):
        """"""
        return self._solver_(*args, **kwargs)

    @classmethod
    def ___parallel_solver_path___(cls):
        return "TOOLS.linear_algebra.solvers.parallel."

    @classmethod
    def ___coded_parallel_solvers___(cls):
        """"""
        return {'GMRES', 'LGMRES', 'BiCGSTAB'}





if __name__ == '__main__':
    solver = ParallelSolverDistributor('BiCGSTAB')
    print(solver._solver_)