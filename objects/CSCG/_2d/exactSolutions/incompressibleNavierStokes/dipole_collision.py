# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/08/29 2:37 PM
See [H.J.H. Clercx, C.-H. Bruneau, The normal and oblique collision of a dipole with a no-slip
boundary, Computers & Fluids, 2006]

and

[G.G. de Diego, A. Palha, M. Gerritsma, Inclusion of no-slip boundary conditions in the MEEVC
scheme, JCP, 2019]

"""
import sys
from abc import ABC

if './' not in sys.path:
    sys.path.append('./')
import numpy as np

from objects.CSCG._2d.exactSolutions.incompressibleNavierStokes.base import \
    incompressibleNavierStokesBase
from objects.CSCG._2d.fields.vector.main import _2dCSCG_VectorField
from components.quadrature import Quadrature
from math import isclose


class DipoleCollision(incompressibleNavierStokesBase, ABC):
    """The 2d dipole collision case. See [H.J.H. Clercx et al., Computers & Fluids, 2006]"""
    def __init__(
            self, mesh, nu=1/625, x1=0, y1=0.1, x2=0, y2=-0.1, r0=0.1, we1=320, we2=-320,
            K_scale=2,
    ):
        """
        The domain must be [-1,1]^2.

        (x1, y1) : center of the monopole 1.
        (x2, y2) : center of the monopole 2.

        These default inputs imply a typical normal collision.

        The velocity is scaled such that initial kinetic energy is equal to 2. And in this case, the initial
        enstrophy will be approximately 800.

        Parameters
        ----------
        mesh
        nu
        x1
        y1
        x2
        y2
        r0
        we1
        we2
        K_scale
        """
        self._x1_ = x1
        self._y1_ = y1
        self._x2_ = x2
        self._y2_ = y2
        self._r0_ = r0
        self._r0_square_ = r0 ** 2
        self._we1_ = we1
        self._we2_ = we2
        assert we1 + we2 == 0
        self._abs_we_ = abs(we1)

        quad_nodes, quad_weights = Quadrature([99, 99], category='Gauss').quad
        quad_nodes = np.meshgrid(*quad_nodes, indexing='ij')

        U0 = self._u_unscaled(*quad_nodes)
        V0 = self._v_unscaled(*quad_nodes)
        K_unscaled = 0.5 * np.einsum('ij, i, j ->', U0**2 + V0**2, quad_weights[0], quad_weights[1], optimize='greedy')

        scale = K_scale / K_unscaled
        self._v_scale_ = np.sqrt(scale)  # 0.9360262042975791, can use approximately 0.936026
        # self._v_scale_ = 0.93603

        super(DipoleCollision, self).__init__(mesh, nu)

        W0 = self.omega(0, *quad_nodes)
        E_scaled = 0.5 * np.einsum('ij, i, j ->', W0**2, quad_weights[0], quad_weights[1], optimize='greedy')
        assert isclose(E_scaled, 800), f"a safety check."

        P0 = self.___palinstrophy_distribution_distribution___(0, *quad_nodes)
        P_scaled = np.einsum('ij, i, j ->', P0, quad_weights[0], quad_weights[1], optimize='greedy')
        assert isclose(441855, P_scaled, abs_tol=1), f"a safety check."
        # print(E_scaled, P_scaled)

        boundary_outline_data = self.mesh.boundaries.visualize.matplot.___outline___(data_only=True)

        assert not mesh.domain.whether.periodic, "DipoleCollision is not for periodic domain."
        for bn in boundary_outline_data:
            for outline in boundary_outline_data[bn]:
                self.___domain_outline_requirements___(*outline)

    @staticmethod
    def ___domain_outline_requirements___(x, y):
        """"""
        U = all(np.logical_and(x == -1, -1 <= y, y <= 1))
        D = all(np.logical_and(x == 1, -1 <= y, y <= 1))
        L = all(np.logical_and(y == -1, -1 <= x, x <= 1))
        R = all(np.logical_and(y == 1, -1 <= x, x <= 1))
        assert U or D or L or R, f"domain needs to be [-1, 1]^2."

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

    @property
    def body_force(self):
        """body force valid all the time."""
        if self._body_force_ is None:
            self._body_force_ = _2dCSCG_VectorField(
                self.mesh, [self.fx, self.fy], valid_time=None, name='body_force')
        return self._body_force_

    def _u_unscaled(self, x, y):
        """"""
        x1 = self._x1_
        y1 = self._y1_
        x2 = self._x2_
        y2 = self._y2_

        r1_square = (x - x1) ** 2 + (y - y1) ** 2
        r2_square = (x - x2) ** 2 + (y - y2) ** 2

        term1 = -0.5 * self._abs_we_ * (y - y1) * np.exp(- r1_square / self._r0_square_)
        term2 = 0.5 * self._abs_we_ * (y - y2) * np.exp(- r2_square / self._r0_square_)

        return term1 + term2

    def u(self, t, x, y):
        return self._v_scale_ * self._u_unscaled(x, y)

    def _v_unscaled(self, x, y):
        x1 = self._x1_
        y1 = self._y1_
        x2 = self._x2_
        y2 = self._y2_

        r1_square = (x - x1) ** 2 + (y - y1) ** 2
        r2_square = (x - x2) ** 2 + (y - y2) ** 2

        term1 = 0.5 * self._abs_we_ * (x - x1) * np.exp(- r1_square / self._r0_square_)
        term2 = -0.5 * self._abs_we_ * (x - x2) * np.exp(- r2_square / self._r0_square_)

        return term1 + term2

    def v(self, t, x, y):
        return self._v_scale_ * self._v_unscaled(x, y)

    def _omega_unscaled(self, x, y):
        """"""
        x1 = self._x1_
        y1 = self._y1_
        x2 = self._x2_
        y2 = self._y2_
        we1 = self._we1_
        we2 = self._we2_

        r1_square = (x - x1) ** 2 + (y - y1) ** 2
        r2_square = (x - x2) ** 2 + (y - y2) ** 2

        omega_1 = we1 * (1 - r1_square / self._r0_square_) * np.exp(- r1_square / self._r0_square_)
        omega_2 = we2 * (1 - r2_square / self._r0_square_) * np.exp(- r2_square / self._r0_square_)

        return omega_2 + omega_1

    def omega(self, t, x, y):
        return self._v_scale_ * self._omega_unscaled(x, y)


if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_2d/exactSolutions/incompressibleNavierStokes/dipole_collision.py
    from objects.CSCG._2d.master import MeshGenerator, ExactSolutionSelector
    mesh = MeshGenerator('crazy', bounds=[[-1, 1], [-1, 1]], c=0.0)([2, 2])
    es = ExactSolutionSelector(mesh)("icpsNS:dipole_collision", show_info=True)

    es.vorticity.visualize(time=0, colormap='bwr')
