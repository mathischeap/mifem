# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('/')
from root.config.main import *
from objects.CSCG._3d.forms.standard.base.main import _3dCSCG_Standard_Form
from objects.CSCG._3d.forms.standard._0_form.special.main import _0Form_Special
from objects.CSCG._3d.forms.standard._0_form.discretize.main import _3dCSCG_Discretize
from objects.CSCG._3d.forms.standard._0_form.reconstruct import _3dCSCG_SF0_reconstruct
from objects.CSCG._3d.forms.standard._0_form.inheriting.private import _3dCSCG_S0F_Private



class _3dCSCG_0Form(_3dCSCG_S0F_Private, _3dCSCG_Standard_Form):
    """
    Standard 0-form.

    :param mesh:
    :param space:
    :param is_hybrid:
    :param orientation:
    :param numbering_parameters:
    :param name:
    """
    def __init__(self, mesh, space, is_hybrid=True,
        orientation='outer', numbering_parameters='Naive',  name=None):
        if name is None:
            if is_hybrid:
                name = 'hybrid-' + orientation + '-oriented-0-form'
            else:
                name = orientation + '-oriented-0-form'
        super().__init__(mesh, space, is_hybrid, orientation, numbering_parameters, name)
        self._k_ = 0
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_standard_0form')
        self._special_ = _0Form_Special(self)
        self.___PRIVATE_reset_cache___()
        self._discretize_ = _3dCSCG_Discretize(self)
        self._reconstruct_ = None
        self._freeze_self_()

    def ___PRIVATE_TW_FUNC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_ScalarField':
            assert func_body.ftype in ('standard',), \
                f"3dCSCG 0form FUNC do not accept func _3dCSCG_ScalarField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 0form FUNC do not accept func {func_body.__class__}")

    def ___PRIVATE_TW_BC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_ScalarField':
            assert func_body.ftype in ('standard','boundary-wise'), \
                f"3dCSCG 0form BC do not accept func _3dCSCG_ScalarField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 0form BC do not accept func {func_body.__class__}")


    def ___PRIVATE_reset_cache___(self):
        super().___PRIVATE_reset_cache___()

    @property
    def special(self):
        return self._special_

    @property
    def discretize(self):
        return self._discretize_

    @property
    def reconstruct(self):
        if self._reconstruct_ is None:
            self._reconstruct_ = _3dCSCG_SF0_reconstruct(self)
        return self._reconstruct_











if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\forms\standard\_0_form\main.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.0)([5,5,5])
    space = SpaceInvoker('polynomials')([('Lobatto',3), ('Lobatto',3), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    # es = ExactSolutionSelector(mesh)('icpsNS:sincosRD')

    def p(t, x, y, z):
        return - 6 * np.pi * np.sin(2*np.pi*x) * np.sin(2*np.pi*y) * np.sin(2*np.pi*z) + 0 * t


    scalar = FC('scalar', p)

    f0 = FC('0-f', is_hybrid=True)
    f0.TW.func.do.set_func_body_as(scalar)
    f0.TW.current_time = 0
    f0.TW.do.push_all_to_instant()
    f0.do.discretize()
    print(f0.error.L())

    df0 = FC('0-adf')
    df0.prime.TW.func.do.set_func_body_as(scalar)
    df0.prime.TW.current_time = 0
    df0.prime.TW.do.push_all_to_instant()
    df0.prime.do.discretize()
    print(df0.prime.error.L())

    # print(f0.___parameters___)

    # regions = mesh.domain.regions['R:R']
    # RS = regions.sub_geometry.make_a_perpendicular_slice_object_on(t=4.5/9)
    # MS = mesh.sub_geometry.make_a_perpendicular_slice_object_on(RS)
    # f0.visualize.matplot.perpendicular_slice(MS, saveto='')

    # f0.export.field.to_file('t0f.mat')
    # print(f0.matrices.trace)