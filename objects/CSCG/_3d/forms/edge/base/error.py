# -*- coding: utf-8 -*-
"""

"""
import sys
if './' not in sys.path: sys.path.append('/')


from components.freeze.main import FrozenOnly
from root.config.main import *




class _3dCSCG_Edge_Error(FrozenOnly):
    """"""
    def __init__(self, ef):
        """"""
        self._ef_ = ef
        self._freeze_self_()


    def L(self, n='infinity', quad_degree=None, quad_density=None):
        """

        :param n:
        :param quad_degree:
        :param quad_density:
        :return:
        """

        assert self._ef_.CF.ftype == 'standard', f"Currently, this L^n error method only works for standard functions."
        assert self._ef_.cochain.local is not None, " I have no cochain."

        quad_degree = [self._ef_.dqp[i] + 2 for i in range(3)] \
            if quad_degree is None else quad_degree


        if n == 'infinity':
            if quad_density is not None:
                assert isinstance(quad_density, (int, float)) and quad_density > 0, \
                    f"quad_density ={quad_density} must be int or float > 0."
                NUM_elements = self._ef_.mesh.elements.global_num
                density_per_element = quad_density / NUM_elements
                num_nodes = density_per_element**(1/3)
                if num_nodes < 1: num_nodes = 3

                if num_nodes % 1 >= 0.5:
                    num_nodes = int(num_nodes) + 1
                else:
                    num_nodes = int(num_nodes)

                _nodes_ = np.linspace(-1, 1, num_nodes+1)
                _nodes_ = (_nodes_[:-1] + _nodes_[1:]) / 2
                quad_nodes = [_nodes_, _nodes_, _nodes_]

            else:
                quad_nodes = [np.linspace(-1, 1, p+10) for p in quad_degree]

        else:
            assert isinstance(n, int) and n > 0, f"L^{n} error is not valid."
            quad_nodes, _, quad_weights = self._ef_.space.___PRIVATE_do_evaluate_quadrature___(quad_degree)

        #-- reconstruct the edge form on all edge-elements -----------------------------------------
        xyz, v = self._ef_.reconstruct(*quad_nodes)

        if n == 'infinity':
            localError = -1
            for i in self._ef_.mesh.edge.elements:

                error_i = [np.max(np.abs(v[i][m] -
                                         self._ef_.CF.___DO_evaluate_func_at_time___()[m](*xyz[i])))
                           for m in range(1)]
                error_i = max(error_i)

                localError = error_i if error_i > localError else localError

            LOC_ERR = COMM.gather(localError, root=MASTER_RANK)

            if RANK == MASTER_RANK:
                globalError = max(LOC_ERR)
            else:
                globalError = None

            globalError = COMM.bcast(globalError, root=MASTER_RANK)
        else:
            raise NotImplementedError(f"Not implemented for L^{n}-error of edge forms.")

        return globalError



if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\forms\edge\base\error.py

    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.)([5,6,7])
    space = SpaceInvoker('polynomials')([('Lobatto',5), ('Lobatto',5), ('Lobatto',5)])
    FC = FormCaller(mesh, space)

    e0 = FC('0-e')

    def p(t, x, y, z): return - 6 * np.pi * np.sin(2*np.pi*x) * np.sin(2*np.pi*y) * np.sin(2*np.pi*z) + 0 * t
    scalar = FC('scalar', p)

    e0.TW.func.do.set_func_body_as(scalar)
    e0.TW.current_time = 0
    e0.TW.do.push_all_to_instant()

    e0.discretize()

    error = e0.error.L(quad_density=100000)
    print(error)