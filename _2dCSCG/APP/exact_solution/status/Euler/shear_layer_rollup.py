
import sys


if './' not in sys.path: sys.path.append('./')

from _2dCSCG.APP.exact_solution.status.Euler.base import EulerBase
from numpy import tanh, sin, pi#, cosh
from screws.miscellaneous.generalized_piecewise_function import genpiecewise

class ShearLayerRollup(EulerBase):
    """
    See Section 5.2 of MEEVC paper.

    """
    def __init__(self, es, delta=pi/15, epsilon=0.05):
        self._delta_ = delta
        self._epsilon_ = epsilon
        super(ShearLayerRollup, self).__init__(es)
        #-------- check the domain ---------------------------------------------------------------
        assert self.mesh.domain.name == 'CrazyPeriodic', \
            f"Shear-Layer-Rollup exact solution only works in crazy_periodic domain, " \
            f"now it is {self.mesh.domain.name}."
        bx, by = self.mesh.domain.domain_input.bounds
        assert tuple(bx) == (0, 2*pi) and tuple(by) == (0, 2*pi), \
            f"ShearLayerRollup can only work in [0, 2pi]^2 periodic domain"
        #-----------------------------------------------------------------------------------------

    @property
    def valid_time(self):
        return 0 # this exact solution will only be valid at t=0

    @property
    def epsilon(self):
        return self._epsilon_

    @property
    def delta(self):
        return self._delta_

    # noinspection PyUnusedLocal
    def _u_low_(self, t, x, y): # y <= pi
        return tanh((y-pi/2)/self.delta)

    # noinspection PyUnusedLocal
    def _u_up_(self, t, x, y): # y > pi
        return tanh((3*pi/2 - y)/self.delta)


    def u(self, t, x, y):
        return genpiecewise([t, x, y], [y <= pi, y > pi], [self._u_low_, self._u_up_])
    def v(self, t, x, y):
        return self.epsilon * sin(x)


    def fx(self, t, x, y):
        return 0 * x
    def fy(self, t, x, y):
        return 0 * y


if __name__ == '__main__':
    # mpiexec -n 4 python _2dCSCG\APP\exact_solution\status\Euler\shear_layer_rollup.py
    from _2dCSCG.main import MeshGenerator, ExactSolutionSelector
    mesh = MeshGenerator('crazy_periodic', bounds=[[0, 2*pi], [0, 2*pi]], c=0.3)([2, 2])
    es = ExactSolutionSelector(mesh)("Euler:shear_layer_rollup", show_info=True)

    es.status.vorticity.visualize(time=0)