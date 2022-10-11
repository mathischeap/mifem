# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 5:43 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from root.config.main import rAnk, mAster_rank, np, cOmm


class miUs_Triangular_SF_Error(FrozenOnly):
    """"""

    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._freeze_self_()

    def L(self, n=2, quad_degree=None):
        """
        The global :math:`L^2` error; it is global, so slaves first send info to the secretary who computes the
        global error and sends it back to all slaves.

        :param int n: (`default`:``2``) :math:`L^{n}` error.
        :param quad_degree: The quadrature degree used to compute the error.
        :return: The global :math:`L^{n}` error.
        :rtype: float
        """
        assert self._sf_.cochain.local is not None, " I have no cochain."
        OneOrTwo = 1 if self._sf_.k in (0, 2) else 2
        quad_degree = [self._sf_.space.p + 2 for _ in range(2)] \
            if quad_degree is None else quad_degree
        quad_nodes, _, quad_weights = self._sf_.space.evaluation.quadrature(quad_degree)
        xi, eta = np.meshgrid(*quad_nodes, indexing='ij')
        xi = xi.ravel('F')
        eta = eta.ravel('F')
        xyz, v = self._sf_.reconstruct(*quad_nodes, ravel=True)

        FUNC = self._sf_.CF.do.evaluate_func_at_time()

        localError = list()

        for i in self._sf_.mesh.elements.indices:
            element = self._sf_.mesh.elements[i]
            detJ = element.coordinate_transformation.Jacobian(xi, eta)
            LEIntermediate = np.sum(
                [(v[i][m] - FUNC[m](*xyz[i]))**n for m in range(OneOrTwo)], axis=0
            )

            localError.append(np.sum(LEIntermediate * detJ * quad_weights))

        core_local = np.sum(localError)
        core_local = cOmm.gather(core_local, root=mAster_rank)

        if rAnk == mAster_rank:
            globalError = np.sum(core_local) ** (1 / n)
        else:
            globalError = None
        globalError = cOmm.bcast(globalError, root=mAster_rank)

        return globalError

    def H(self, d_func=None, quad_degree=None):
        """
        Global :math:`H^1`-error; :math:`H^1` error includes :math:`H(\mathrm{curl})` and :math:`H(\mathrm{div})`;
        since it basically use the ``globalL2`` method, so the computation is done in the secretary core and
        communications are needed.

        :param d_func: The function of the derivative of ``self`` form.
        :param quad_degree: The quadrature degree used to compute the error.
        :return: The global :math:`H^{1}` error.
        :rtype: float
        """
        selfErrorL2 = self.L(n=2, quad_degree=quad_degree)
        D_self = self._sf_.coboundary()
        if d_func is None:
            func = self._sf_.CF
            if self._sf_.__class__.__name__ == 'miUsTriangular_S0F_Inner':
                d_func = func.numerical.grad
            elif self._sf_.__class__.__name__ == 'miUsTriangular_S0F_Outer':
                d_func = func.numerical.curl
            elif self._sf_.__class__.__name__ == 'miUsTriangular_S1F_Inner':
                d_func = func.numerical.rot
            elif self._sf_.__class__.__name__ == 'miUsTriangular_S1F_Outer':
                d_func = func.numerical.div
            else:
                raise Exception(f"{self._sf_.__class__.__name__} has no coboundary.")
        else:
            pass

        D_self.CF = d_func
        D_self.CF.current_time = self._sf_.CF.current_time
        DErrorL2 = D_self.error.L(n=2, quad_degree=quad_degree)
        return (selfErrorL2 ** 2 + DErrorL2 ** 2) ** 0.5


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
