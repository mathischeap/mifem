# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/08/29 7:29 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from tools.linearAlgebra.elementwiseCache.objects.columnVector.main import EWC_ColumnVector
from objects.CSCG._2d.forms.standard._0_form.outer.boundary_integrate.helpers.ipS_helper import ipS_Helper


class _2dCSCG_Outer_S0Form_BI(FrozenOnly):
    """"""

    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._freeze_self_()


    def inner_product_with(self, S, quad_degree=None):
        """Let s1f be denoted as w. We do (w,  S)_{\partial\Omega} here. S must be a given scalar.
        And we get a vector.

        This returns a `EWC_ColumnVector` whose local vector refers to the
        local dofs (basis functions) of the self form, i.e., `s1f`.

        Parameters
        ----------
        S : must be a Scalar
        quad_degree :

        Returns
        -------
        EWCv : EWC_ColumnVector

        """

        assert S.mesh == self._sf_.mesh, f"meshes do not match."

        assert S.__class__.__name__ == '_2dCSCG_ScalarField', f"I need a scalar."

        VDG = ipS_Helper(self._sf_, S, quad_degree)

        # no cache, vector.current_time may change ----------
        return EWC_ColumnVector(self._sf_.mesh, VDG, 'no_cache')


if __name__ == "__main__":
    # mpiexec -n 4 python objects/CSCG/_2d/forms/standard/_0_form/outer/boundary_integrate/main.py
    from objects.CSCG._2d.master import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.)([2, 2])
    space = SpaceInvoker('polynomials')([1,1])
    FC = FormCaller(mesh, space)

    f = FC('0-f-o', is_hybrid=False)

    BV = {'Upper': 0,
          'Down': 0,
          'Left': 0,
          'Right': 1,}

    S = FC('scalar', BV, name='boundary-vector')
    S.current_time = 0

    B = f.do.boundary_integrate.inner_product_with(S)


    for i in B:
        print(i, B._DG_(i))