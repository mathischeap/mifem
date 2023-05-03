# -*- coding: utf-8 -*-
"""
Yi Zhang
zhangyi_aero@hotmail.com
created at: 3/17/2023 1:45 PM
"""
import sys
from abc import ABC

if './' not in sys.path:
    sys.path.append('./')

from objects.CSCG._2d.exactSolutions.incompressibleNavierStokes.base import \
    incompressibleNavierStokesBase
from objects.CSCG._2d.fields.vector.main import _2dCSCG_VectorField
from components.miscellaneous.generalized_piecewise_function import genpiecewise


class FlowRoundCylinder(incompressibleNavierStokesBase, ABC):
    """"""
    def __init__(
            self, mesh, nu=1/1000, inlet_velocity=1, outlet_pressure=0,
    ):
        """"""
        self._inlet_velocity = inlet_velocity
        self._outlet_pressure = outlet_pressure
        di = mesh.domain._domain_input_
        self._velocity_bc_ = None
        assert di.__class__.__name__ == 'CylinderInChannel'
        self._r = di._r_
        self._li = di._li_
        super(FlowRoundCylinder, self).__init__(mesh, nu)

    @property
    def valid_time(self):
        return "valid_only_at_its_first_instant"  # usually t = 0.

    @property
    def Re(self):
        return 1 / self.nu

    def fx(self, t, x, y):
        return 0 * x

    def fy(self, t, x, y):
        return 0 * y


    def u(self, t, x, y):
        return 0 * x

    def v(self, t, x, y):
        return 0 * x

    def _u_bc_0(self, t, x, y):
        return 0 * x + self._inlet_velocity

    def _u_bc_1(self, t, x, y):
        return 0 * x

    def _v_bc(self, t, x, y):
        return 0 * x

    def _u_bc(self, t, x, y):
        return genpiecewise([t, x, y], [x < -(self._r+self._li), x >= -(self._r+self._li)], [self._u_bc_0, self._u_bc_1])

    @property
    def velocity_bc(self):
        """A scalar field of the kinetic energy distribution."""
        if self._velocity_bc_ is None:
            self._velocity_bc_ = _2dCSCG_VectorField(
                self.mesh, [self._u_bc, self._v_bc], valid_time=self.valid_time, name='velocity')
        return self._velocity_bc_


    def _tp_(self, t, x, y):  # total pressure
        return 0 * x



if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_2d/exactSolutions/incompressibleNavierStokes/flow_round_cylinder.py
    from objects.CSCG._2d.master import MeshGenerator, ExactSolutionSelector
    mesh = MeshGenerator('cic')([2, 2])
    es = ExactSolutionSelector(mesh)("icpsNS:flow_round_cylinder", show_info=True)

    es.velocity.visualize(time=0, colormap='bwr')