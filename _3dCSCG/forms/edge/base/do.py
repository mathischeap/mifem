
"""

"""
import sys
if './' not in sys.path: sys.path.append('../')


from screws.frozen import FrozenOnly
from root.config import *




class _3dCSCG_Edge_DO(FrozenOnly):
    """"""
    def __init__(self, ef):
        """"""
        self._ef_ = ef
        self._freeze_self_()



    def evaluate_basis_at_meshgrid(self, xi, eta, sigma):
        return self._ef_.space.do.evaluate_edge_basis_at_meshgrid(self._ef_.k, xi, eta, sigma)






if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\form\edge\do.py
    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.25)([5,6,7])
    space = SpaceInvoker('polynomials')([('Lobatto',2), ('Lobatto',1), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    e0 = FC('0-e')

    def p(t, x, y, z): return - 6 * np.pi * np.sin(2*np.pi*x) * np.sin(2*np.pi*y) * np.sin(2*np.pi*z) + 0 * t
    scalar = FC('scalar', p)

    e0.TW.func.do.set_func_body_as(scalar)
    e0.TW.current_time = 0
    e0.TW.do.push_all_to_instant()

    e0.discretize()