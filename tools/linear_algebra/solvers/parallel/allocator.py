"""
A collection of all parallel solvers.
"""



import sys
if './' not in sys.path: sys.path.append('./')


from importlib import import_module


from screws.freeze.main import FrozenOnly


class ParallelSolverDistributor(FrozenOnly):
    """"""
    def __init__(self, solver_ID, routine='auto', name="no-name-solver"):
        """"""
        solver_name = self.___coded_parallel_solvers___()[solver_ID]
        solver_path = self.___parallel_solver_path___()[solver_ID ]
        self._solver_ = getattr(import_module(solver_path), solver_name)(routine, name)

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
        assert A.__class__.__name__ == "GlobalMatrix", \
            f"A needs to be a GlobalMatrix'. Now I get {A.__class__}."
        assert b.__class__.__name__ == "GlobalVector", \
            f"b needs to be a 'GlobalVector'. Now I get {b.__class__}."
        Results = self._solver_(A, b, *args, **kwargs)
        assert Results[0].__class__.__name__ == 'LocallyFullVector', f"results must be a LocallyFullVector!"
        return Results
    @classmethod
    def ___coded_parallel_solvers___(cls):
        """"""
        return {'GMRES'   : 'GMRES',
                'LGMRES'  : 'LGMRES',
                'BiCGSTAB': 'BiCGSTAB',
                'direct'  : 'Direct'}


    @classmethod
    def ___parallel_solver_path___(cls):
        base_path = '.'.join(str(cls).split(' ')[1][1:-2].split('.')[:-2]) + '.'
        return {
            'GMRES'   : base_path + 'GMRES.main',
            'LGMRES'  : base_path + 'LGMRES.main',
            'BiCGSTAB': base_path + 'BiCGSTAB.main',
            'direct'  : base_path + 'direct.main'
        }





if __name__ == '__main__':
    solver = ParallelSolverDistributor('BiCGSTAB')
    print(solver._solver_)