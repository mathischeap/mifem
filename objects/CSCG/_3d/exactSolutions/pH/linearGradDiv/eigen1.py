# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/23/2022 5:37 PM
"""
import sys
from numpy import sqrt, sin, cos

if './' not in sys.path: sys.path.append('./')
from objects.CSCG._3d.exactSolutions.pH.linearGradDiv.base import pH_LinearGradDiv_Base



class Eigen1(pH_LinearGradDiv_Base):
    """The first eigen-solution.
    """

    def __init__(self, mesh, OM=(1, 1, 1), PHI=(0, 0, 0, 0)):
        """"""
        self._OM_ = OM    # om_x, om_y, om_z
        self._PHI_ = PHI  # phi_x, phi_y, phi_z, phi_t
        super(Eigen1, self).__init__(mesh) # freeze and check here!

    def p(self, t, x, y, z):
        """Translated from Andrea's Firedrake codes.

        Parameters
        ----------
        t
        x
        y
        z

        Returns
        -------

        """
        om_x, om_y, om_z = self._OM_
        om_t = sqrt(om_x ** 2 + om_y ** 2 + om_z ** 2)
        phi_x, phi_y, phi_z, phi_t = self._PHI_

        dft = om_t * (2 * cos(om_t * t + phi_t) - 3 * sin(om_t * t + phi_t))
        g_xyz = cos(om_x * x + phi_x) * sin(om_y * y + phi_y) * sin(om_z * z + phi_z)

        return g_xyz * dft

    def u(self, t, x, y, z):
        om_x, om_y, om_z = self._OM_
        om_t = sqrt(om_x ** 2 + om_y ** 2 + om_z ** 2)
        phi_x, phi_y, phi_z, phi_t = self._PHI_

        ft = 2 * sin(om_t * t + phi_t) + 3 * cos(om_t * t + phi_t)
        dg_xyz_x = - om_x * sin(om_x * x + phi_x) * sin(om_y * y + phi_y) * sin(om_z * z + phi_z)

        return dg_xyz_x * ft

    def v(self, t, x, y, z):
        om_x, om_y, om_z = self._OM_
        om_t = sqrt(om_x ** 2 + om_y ** 2 + om_z ** 2)
        phi_x, phi_y, phi_z, phi_t = self._PHI_

        ft = 2 * sin(om_t * t + phi_t) + 3 * cos(om_t * t + phi_t)
        dg_xyz_y = om_y * cos(om_x * x + phi_x) * cos(om_y * y + phi_y) * sin(om_z * z + phi_z)

        return dg_xyz_y * ft

    def w(self, t, x, y, z):
        om_x, om_y, om_z = self._OM_
        om_t = sqrt(om_x ** 2 + om_y ** 2 + om_z ** 2)
        phi_x, phi_y, phi_z, phi_t = self._PHI_

        ft = 2 * sin(om_t * t + phi_t) + 3 * cos(om_t * t + phi_t)
        dg_xyz_z = om_z * cos(om_x * x + phi_x) * sin(om_y * y + phi_y) * cos(om_z * z + phi_z)

        return dg_xyz_z * ft





if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_3d/exactSolutions/pH/linearGradDiv/eigen1.py
    from objects.CSCG._3d.master import MeshGenerator

    mesh = MeshGenerator('crazy', c=0.0)([5, 5, 5])

    es = Eigen1(mesh)
