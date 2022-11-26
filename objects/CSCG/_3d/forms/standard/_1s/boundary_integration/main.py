# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/26/2022 9:38 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from components.freeze.main import FrozenOnly

from tools.linearAlgebra.elementwiseCache.objects.columnVector.main import EWC_ColumnVector
from objects.CSCG._3d.forms.standard._1s.boundary_integration.helpers.V_helper import S1F_BI_V_Helper


class _3dCSCG_S1F_BI(FrozenOnly):
    """"""

    def __init__(self, s1f):
        """"""
        self._s1f_ = s1f
        self._freeze_self_()

    def inner_product_with(self, V, quad_degree=None):
        """Let s1f be denoted as w. We do (w,  V)_{\partial\Omega} here. V must be given. And we
        get a vector.

        This returns a `EWC_ColumnVector` whose local vector refers to the
        local dofs (basis functions) of the self form, i.e., `s1f`.

        Parameters
        ----------
        V : must be a vector
        quad_degree :

        Returns
        -------
        EWCv : EWC_ColumnVector

        """

        assert V.mesh == self._s1f_.mesh, f"meshes do not match."

        assert V.__class__.__name__ == '_3dCSCG_VectorField', f"I need a vector."

        VDG = S1F_BI_V_Helper(self._s1f_, V, quad_degree)

        # no cache, vector.current_time may change
        return EWC_ColumnVector(self._s1f_.mesh, VDG, 'no_cache')








if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_3d/forms/standard/_1s/boundary_integration/main.py

    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.)([10,10,10])
    space = SpaceInvoker('polynomials')([1,1,1])
    FC = FormCaller(mesh, space)

    f1 = FC('1-f', is_hybrid=False)

    # noinspection PyUnusedLocal
    def ZERO(t, x, y, z): return 0 * x
    # noinspection PyUnusedLocal
    def ONE(t, x, y, z): return 1 + 0 * x

    BV = {'North': [ZERO, ZERO, ZERO],
          'South': [ZERO, ZERO, ZERO],
          'West': [ZERO, ZERO, ZERO],
          'East': [ZERO, ZERO, ZERO],
          'Back': [ZERO, ZERO, ZERO],
          'Front': [ONE, ZERO, ZERO],}

    V = FC('vector', BV, name='boundary-vector')
    VP = V.components.T_perp
    VP.current_time = 0

    B = f1.do.boundary_integrate.inner_product_with(VP)


    for i in B:
        print(i, B[i])