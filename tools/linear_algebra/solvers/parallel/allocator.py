"""
A collection of all parallel solvers.
"""



import sys
if './' not in sys.path: sys.path.append('./')


from importlib import import_module


from screws.freeze.main import FrozenOnly


class ParallelSolverDistributor(FrozenOnly):
    """"""
    def __init__(self, solver_name, routine='auto', name="no-name-solver"):
        """"""
        self._sn_ = solver_name
        solver_path = self.___parallel_solver_path___() + solver_name + '.main'
        self._solver_ = getattr(import_module(solver_path), solver_name)(routine, name)

    def __call__(self, A, b,  *args, **kwargs):
        """"""
        assert A.__class__.__name__ == "GlobalMatrix", \
            f"A needs to be a GlobalMatrix'. Now I get {A.__class__}."
        assert b.__class__.__name__ == "GlobalVector", \
            f"b needs to be a 'GlobalVector'. Now I get {b.__class__}."
        return self._solver_(A, b, *args, **kwargs)

    @classmethod
    def ___parallel_solver_path___(cls):
        return "tools.linear_algebra.solvers.parallel."

    @classmethod
    def ___coded_parallel_solvers___(cls):
        """"""
        return {'GMRES', 'LGMRES', 'BiCGSTAB'}





if __name__ == '__main__':
    solver = ParallelSolverDistributor('BiCGSTAB')
    print(solver._solver_)