# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('/')
from objects.CSCG._3d.forms.standard._3s.discretize.main import _3dCSCG_Discretize
from objects.CSCG._3d.forms.standard.base.main import _3dCSCG_Standard_Form
from objects.CSCG._3d.forms.standard._3s.special import _3Form_Special
from objects.CSCG._3d.forms.standard._3s.reconstruct import _3dCSCG_SF3_Reconstruct
from objects.CSCG._3d.forms.standard._3s.inheriting.private import _3dCSCG_S3F_Private



class _3dCSCG_3Form(_3dCSCG_S3F_Private, _3dCSCG_Standard_Form):
    """Standard 3-form.

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
                name = 'hybrid-' + orientation + '-oriented-3-form'
            else:
                name = orientation + '-oriented-3-form'
        super().__init__(mesh, space, is_hybrid, orientation, numbering_parameters, name)
        self._k_ = 3
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_standard_3form')
        self._special_ = _3Form_Special(self)
        self.___PRIVATE_reset_cache___()
        self._discretize_ = _3dCSCG_Discretize(self)
        self._reconstruct_ = None
        self._freeze_self_()

    def ___PRIVATE_reset_cache___(self):
        super().___PRIVATE_reset_cache___()

    def ___PRIVATE_TW_FUNC_body_checker___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_ScalarField':
            assert func_body.ftype in ('standard',), \
                f"3dCSCG 3form FUNC do not accept func _3dCSCG_ScalarField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 3form FUNC do not accept func {func_body.__class__}")

    @property
    def special(self):
        return self._special_

    @property
    def discretize(self):
        return self._discretize_

    @property
    def reconstruct(self):
        if self._reconstruct_ is None:
            self._reconstruct_ = _3dCSCG_SF3_Reconstruct(self)
        return self._reconstruct_





if __name__ == '__main__':
    # mpiexec python _3dCSCG\forms\standard\_3_form\main.py

    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.25)([3,3,3])
    space = SpaceInvoker('polynomials')([('Lobatto',2), ('Lobatto',2), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    es = ExactSolutionSelector(mesh)('icpsNS:sincosRD')
    f3 = FC('3-f', is_hybrid=False)

    f3.TW.func.do.set_func_body_as(es, 'pressure')
    f3.TW.current_time = 0
    f3.TW.___DO_push_all_to_instant___()
    f3.do.discretize()

    # from tools.CSCG.partial_dofs import PartialDofs
    #
    # pd = PartialDofs(f3)