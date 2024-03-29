# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/26/2022 4:36 PM
"""
import sys

if './' not in sys.path:
    sys.path.append('./')


from objects.CSCG._3d.exactSolutions.incompressibleNavierStokes.base import \
    incompressible_NavierStokes_Base

from objects.CSCG._3d.fields.vector.main import _3dCSCG_VectorField


# noinspection PyUnusedLocal
def _0_(t, x, y, z): return 0 * x


class LidDrivenCavity(incompressible_NavierStokes_Base):
    """The classic lid driven cavity case."""

    def __init__(self, mesh, nu=0.01, lid_velocity=1):
        """"""
        self._lid_velocity_ = lid_velocity
        super(LidDrivenCavity, self).__init__(mesh, nu)

    # noinspection PyUnusedLocal
    def _V_(self, t, x, y, z):
        return self._lid_velocity_ + 0 * x

    @property
    def lid_velocity(self):
        return self._lid_velocity_

    @property
    def velocity(self):
        """The boundary velocity. We should only use this for the boundary condition.

        And this boundary velocity is valid all time.
        """
        if self._velocity_ is None:
            BV = {'North': [_0_, _0_, _0_],
                  'South': [_0_, _0_, _0_],
                  'West':  [_0_, _0_, _0_],
                  'East':  [_0_, _0_, _0_],
                  'Back':  [_0_, _0_, _0_],
                  'Front': [self._V_, _0_, _0_]}

            self._velocity_ = _3dCSCG_VectorField(
                self.mesh,
                BV,
                ftype='boundary-wise',
                valid_time=None,
                name='velocity_boundary_condition'
                )

        return self._velocity_

    def fx(self, t, x, y, z): return 0 * x
    def fy(self, t, x, y, z): return 0 * x
    def fz(self, t, x, y, z): return 0 * x

    def ___PreFrozenChecker___(self):
        """Skip checking."""
        pass


if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_3d/exact_solutions/status/incompressibleNavierStokes/lid_driven_cavity.py
    from objects.CSCG._3d.master import MeshGenerator, ExactSolutionSelector
    mesh = MeshGenerator('crazy', c=0.0)([2, 2, 2])
    es = ExactSolutionSelector(mesh)("icpsNS:LDC", show_info=True)

    T_perp = es.status.velocity.components.T_perp

    T_perp.current_time = 0
    T_perp.visualize()
