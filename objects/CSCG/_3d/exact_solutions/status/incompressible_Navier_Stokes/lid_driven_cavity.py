# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/26/2022 4:36 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')


from objects.CSCG._3d.exact_solutions.status.incompressible_Navier_Stokes.base import \
    incompressible_NavierStokes_Base

from objects.CSCG._3d.fields.vector.main import _3dCSCG_VectorField



class LidDrivenCavity(incompressible_NavierStokes_Base):
    """"""

    def __init__(self, es, nu=0.01, lid_velocity=1):
        """"""
        self._lid_velocity_ = lid_velocity
        super(LidDrivenCavity, self).__init__(es, nu)

    # noinspection PyUnusedLocal
    def _0_(self, t, x, y, z): return 0 * x

    # noinspection PyUnusedLocal
    def _V_(self, t, x, y, z): return self._lid_velocity_ + 0 * x

    @property
    def lid_velocity(self):
        return self._lid_velocity_

    @property
    def velocity(self):
        if self._velocity_ is None:
            BV = {'North': [self._0_, self._0_, self._0_],
                  'South': [self._0_, self._0_, self._0_],
                  'West':  [self._0_, self._0_, self._0_],
                  'East':  [self._0_, self._0_, self._0_],
                  'Back':  [self._0_, self._0_, self._0_],
                  'Front': [self._V_, self._0_, self._0_]}

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
    # mpiexec -n 4 python objects/CSCG/_3d/exact_solutions/status/incompressible_Navier_Stokes/lid_driven_cavity.py
    from objects.CSCG._3d.master import MeshGenerator, ExactSolutionSelector
    mesh = MeshGenerator('crazy', c=0.0)([2, 2, 2])
    es = ExactSolutionSelector(mesh)("icpsNS:LDC", show_info=True)

    print(es.status.velocity)
    print(es.status.body_force)

