# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/08/29 7:29 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from components.freeze.base import FrozenOnly
from tools.elementwiseCache.dataStructures.objects.columnVector.main import EWC_ColumnVector
from objects.CSCG._2d.forms.standard._1_form.outer.boundary_integrate.helpers.ipVhelper import ipV_Helper


class _2dCSCG_Outer_S1Form_BI(FrozenOnly):
    """"""

    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._freeze_self_()


    def inner_product_with(self, V, quad_degree=None):
        """Let s1f be denoted as w. We do (w,  V)_{\partial\Omega} here. V must be a given vector.
        And we get a vector.

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

        assert V.mesh == self._sf_.mesh, f"meshes do not match."

        assert V.__class__.__name__ == '_2dCSCG_VectorField', f"I need a vector."

        VDG = ipV_Helper(self._sf_, V, quad_degree)

        # no cache, vector.current_time may change
        return EWC_ColumnVector(self._sf_.mesh, VDG, 'no_cache')


if __name__ == "__main__":
    # mpiexec -n 4 python objects/CSCG/_2d/forms/standard/_1_form/outer/boundary_integrate/main.py
    from objects.CSCG._2d.master import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.)([2, 2])
    space = SpaceInvoker('polynomials')([1,1])
    FC = FormCaller(mesh, space)

    f1 = FC('1-f-o', is_hybrid=False)

    BV = {'Upper': [0, 0],
          'Down': [0, 0],
          'Left': [0, 0],
          'Right': [1, 0],}

    V = FC('vector', BV, name='boundary-vector')
    V.current_time = 0

    B = f1.do.boundary_integrate.inner_product_with(V)


