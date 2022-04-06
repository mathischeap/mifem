# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('./')
from objects.CSCG._2d.forms.standard._2_form.inner.special.main import _2Form_Inner_Special
from objects.CSCG._2d.forms.standard._2_form.base.main import _2Form_BASE

class _2dCSCG_2Form_Inner(_2Form_BASE):
    """
    Inner standard 2-form.

    :param mesh:
    :param space:
    :param is_hybrid:
    :param numbering_parameters:
    :param name:
    """
    def __init__(self, mesh, space, is_hybrid=True,
        numbering_parameters='Naive',  name='inner-oriented-2-form'):
        super().__init__(mesh, space, is_hybrid, 'inner', numbering_parameters, name)
        super().__init_2form_base__()
        self._k_ = 2
        self.standard_properties.___PRIVATE_add_tag___('2dCSCG_standard_inner_2form')
        self.standard_properties.___PRIVATE_add_tag___('2dCSCG_standard_2form')
        self._special_ = _2Form_Inner_Special(self)
        self.___PRIVATE_reset_cache___()
        self._freeze_self_()

    @property
    def special(self):
        return self._special_








if __name__ == '__main__':
    # mpiexec python _2dCSCG\form\standard\_2_form_inner.py

    from objects.CSCG._2d.master import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.3)([50,45])
    # mesh = MeshGenerator('chp1',)([2,2])
    space = SpaceInvoker('polynomials')([('Lobatto',3), ('Lobatto',4)])
    FC = FormCaller(mesh, space)

    ES = ExactSolutionSelector(mesh)('sL:sincos1')

    f2 = FC('2-f-i', is_hybrid=False, name='potential')
    f2.TW.func.do.set_func_body_as(ES, 'potential')
    f2.TW.current_time = 0
    f2.TW.do.push_all_to_instant()
    f2.discretize()

    f2.visualize()
    # print(f2.error.L())
    #
    # # from root.mifem import save
    # #
    # # save(f2, 'test_2d_f2_i')
    #
    #
    # from tools.CSCG.partial_cochain.partial_dofs.main import PartialDofs
    #
    # pd = PartialDofs(f2)
