# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/26/2022 2:19 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from objects.CSCG._3d.forms.base import _3dCSCG_FORM_BASE

from objects.CSCG._3d.forms.localTrace.base.cochain.main import _3dCSCG_LocalTrace_Cochain
from objects.CSCG._3d.forms.localTrace.base.numbering.main import _3dCSCG_LocalTrace_Numbering
from objects.CSCG._3d.forms.localTrace.base.num import _3dCSCG_LocalTrace_NUM
from objects.CSCG._3d.forms.localTrace.base.do import _3dCSCG_LocalTrace_Do
from objects.CSCG._3d.forms.localTrace.base.matrices.main import _3dCSCG_LocalTrace_Matrices
from objects.CSCG._3d.forms.localTrace.base.operators.main import _3dCSCG_LocalTrace_Operators
from objects.CSCG._3d.forms.localTrace.base.whether import _3dCSCG_LocalTrace_Whether

class _3dCSCG_LocalTrace(_3dCSCG_FORM_BASE, ndim=3):
    """"""

    def __init__(self, mesh, space, hybrid, orientation, numbering_parameters, name):
        """

        Parameters
        ----------
        mesh
        space
        hybrid : bool
            Whether the local-trace-form is locally in each mesh-element hybrid? Note that this
            does not refer to whether across mesh-elements the form is hybrid.
        orientation
        numbering_parameters
        name
        """
        super().__init__(mesh, space, name)
        assert isinstance(hybrid, bool), f"hybrid be True of False"
        assert orientation in ('inner', 'outer'), " orientation needs to be 'inner' or 'outer'."
        self._orientation_ = orientation
        self._hybrid_ = hybrid
        self._whether_ = _3dCSCG_LocalTrace_Whether(self)
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_local_trace_form')

        self._numbering_ = _3dCSCG_LocalTrace_Numbering(self, numbering_parameters)
        self._cochain_ = _3dCSCG_LocalTrace_Cochain(self)
        self._num_ = _3dCSCG_LocalTrace_NUM(self)
        self._do_ = _3dCSCG_LocalTrace_Do(self)
        self._operators_ = None
        self._matrices_ = None

        self._reconstruct_ = None
        self._discretize_ = None
        self._visualize_ = None

    @property
    def orientation(self):
        return self._orientation_

    @property
    def numbering(self):
        return self._numbering_

    @property
    def cochain(self):
        """Collections of all cochain-related sub-properties."""
        return self._cochain_

    @property
    def num(self):
        return self._num_

    @property
    def do(self):
        return self._do_

    @property
    def operators(self):
        if self._operators_ is None:
            self._operators_ = _3dCSCG_LocalTrace_Operators(self)
        return self._operators_

    @property
    def matrices(self):
        if self._matrices_ is None:
            self._matrices_ = _3dCSCG_LocalTrace_Matrices(self)
        return self._matrices_

    @property
    def whether(self):
        return self._whether_

    @property
    def reconstruct(self):
        return self._reconstruct_

    @property
    def discretize(self):
        return self._discretize_

    @property
    def visualize(self):
        return self._visualize_

if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_3d/forms/localTrace/base/main.py

    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.0)([2,2,2])
    space = SpaceInvoker('polynomials')([3,3,3])
    FC = FormCaller(mesh, space)

    lt0 = FC('0-lt')

    GM0 = lt0.numbering.gathering