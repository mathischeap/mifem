# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
from abc import ABC

if './' not in sys.path:
    sys.path.append('./')

from objects.CSCG._3d.forms.standard._2s.discretize.main import _3dCSCG_Discretize
from objects.CSCG._3d.forms.standard.base.main import _3dCSCG_Standard_Form
from objects.CSCG._3d.forms.standard._2s.special.main import _2Form_Special
from objects.CSCG._3d.forms.standard._2s.project.main import _2Form_Projection
from objects.CSCG._3d.forms.standard._2s.reconstruct import _3dCSCG_SF2_reconstruct
from objects.CSCG._3d.forms.standard._2s.inheriting.private import _3dCSCG_S2F_Private
from objects.CSCG._3d.forms.standard._2s.visualize.main import _3dCSCG_S2F_VISUALIZE
from objects.CSCG._3d.forms.standard._2s.do.main import _3dCSCG_S2F_Do


class _3dCSCG_2Form(_3dCSCG_S2F_Private, _3dCSCG_Standard_Form, ABC):
    """
    Standard 2-form.

    :param mesh:
    :param space:
    :param hybrid:
    :param orientation:
    :param numbering_parameters:
    :param name:
    """
    def __init__(
            self, mesh, space, hybrid=True,
            orientation='outer', numbering_parameters='Naive',  name=None
    ):
        if name is None:
            if hybrid:
                name = 'hybrid-' + orientation + '-oriented-2-form'
            else:
                name = orientation + '-oriented-2-form'
        super().__init__(mesh, space, hybrid, orientation, numbering_parameters, name)
        self._k_ = 2
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_standard_2form')
        self._special_ = _2Form_Special(self)
        self._projection_ = _2Form_Projection(self)
        self._discretize_ = None
        self._reconstruct_ = None
        self._visualize_ = None
        self._DO_ = _3dCSCG_S2F_Do(self)
        self.___kwargs___ = {
            'hybrid': hybrid,
            'orientation': orientation,
            'numbering_parameters': numbering_parameters,
            'name': name,
        }
        self._freeze_self_()

    def ___Pr_check_CF___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_VectorField':
            assert func_body.ftype in ('standard',), \
                f"3dCSCG 2form FUNC do not accept func _3dCSCG_VectorField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 2form FUNC do not accept func {func_body.__class__}")

    def ___Pr_check_BC_CF___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_VectorField':
            assert func_body.ftype in ('standard', 'boundary-wise'), \
                f"3dCSCG 2form BC do not accept func _3dCSCG_VectorField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 2form BC do not accept func {func_body.__class__}")

    @property
    def special(self):
        return self._special_

    @property
    def projection(self):
        """A wrapper of all projection methods."""
        return self._projection_

    @property
    def discretize(self):
        if self._discretize_ is None:
            self._discretize_ = _3dCSCG_Discretize(self)
        return self._discretize_

    @property
    def reconstruct(self):
        if self._reconstruct_ is None:
            self._reconstruct_ = _3dCSCG_SF2_reconstruct(self)
        return self._reconstruct_

    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = _3dCSCG_S2F_VISUALIZE(self)
        return self._visualize_


if __name__ == '__main__':
    # mpiexec -n 5 python objects/CSCG/_3d/forms/standard/_2s/main.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller

    K = 3
    N = 2
    mesh = MeshGenerator('crazy', c=0.0)([K, K, K])
    space = SpaceInvoker('polynomials')([2, 2, 2])
    FC = FormCaller(mesh, space)

    from numpy import cos, sin, pi

    def u(t, x, y, z): return 2 * pi * cos(2*pi*x) * sin(2*pi*y) * sin(2*pi*z) + 0 * t
    def v(t, x, y, z): return 2 * pi * sin(2*pi*x) * cos(2*pi*y) * sin(2*pi*z) + 0 * t
    def w(t, x, y, z): return 2 * pi * sin(2*pi*x) * sin(2*pi*y) * cos(2*pi*z) + 0 * t

    velocity = FC('vector', (u, v, w))

    u = FC('2-f', hybrid=False, name='velocity')

    u.CF = velocity
    u.CF.current_time = 0
    u.discretize()

    xyz, value = u.reconstruct([-1, 0, 1], [-1, 0, 1], [-1, 0, 1])
    if 5 in value:
        U, V, W = value[5]
        print(U.ravel('F'))
        print(V.ravel('F'))
        print(W.ravel('F'))
