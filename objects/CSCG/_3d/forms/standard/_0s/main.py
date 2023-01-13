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
from objects.CSCG._3d.forms.standard.base.main import _3dCSCG_Standard_Form
from objects.CSCG._3d.forms.standard._0s.special.main import _0Form_Special
from objects.CSCG._3d.forms.standard._0s.discretize.main import _3dCSCG_Discretize
from objects.CSCG._3d.forms.standard._0s.reconstruct import _3dCSCG_SF0_reconstruct
from objects.CSCG._3d.forms.standard._0s.inheriting.private import _3dCSCG_S0F_Private
from objects.CSCG._3d.forms.standard._0s.visualize.main import _3dCSCG_S0F_VISUALIZE
from objects.CSCG._3d.forms.standard._0s.do.main import _3dCSCG_S0F_Do


class _3dCSCG_0Form(_3dCSCG_S0F_Private, _3dCSCG_Standard_Form, ABC):
    """
    Standard 0-form.

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
                name = 'hybrid-' + orientation + '-oriented-0-form'
            else:
                name = orientation + '-oriented-0-form'
        super().__init__(mesh, space, hybrid, orientation, numbering_parameters, name)
        self._k_ = 0
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_standard_0form')
        self._special_ = _0Form_Special(self)
        self._discretize_ = None
        self._reconstruct_ = None
        self._visualize_ = None
        self._DO_ = _3dCSCG_S0F_Do(self)
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

        if func_body.__class__.__name__ == '_3dCSCG_ScalarField':
            assert func_body.ftype in ('standard',), \
                f"3dCSCG 0form FUNC do not accept func _3dCSCG_ScalarField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 0form FUNC do not accept func {func_body.__class__}")

    def ___Pr_check_BC_CF___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_ScalarField':
            assert func_body.ftype in ('standard', 'boundary-wise'), \
                f"3dCSCG 0form BC do not accept func _3dCSCG_ScalarField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 0form BC do not accept func {func_body.__class__}")

    @property
    def special(self):
        return self._special_

    @property
    def discretize(self):
        if self._discretize_ is None:
            self._discretize_ = _3dCSCG_Discretize(self)
        return self._discretize_

    @property
    def reconstruct(self):
        if self._reconstruct_ is None:
            self._reconstruct_ = _3dCSCG_SF0_reconstruct(self)
        return self._reconstruct_

    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = _3dCSCG_S0F_VISUALIZE(self)
        return self._visualize_


if __name__ == '__main__':
    # mpiexec -n 6 python objects/CSCG/_3d/forms/standard/_0s/main.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller

    mesh = MeshGenerator('crazy', c=0.)([3, 3, 3])
    space = SpaceInvoker('polynomials')([2, 2, 2])
    FC = FormCaller(mesh, space)

    f0 = FC('3-f', hybrid=True)

    M = f0.matrices.mass

    print(M.condition.pseudo_sparsity)
