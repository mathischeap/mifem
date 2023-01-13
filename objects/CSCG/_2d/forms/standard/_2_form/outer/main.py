# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path:
    sys.path.append('./')
from objects.CSCG._2d.forms.standard._2_form.outer.special import _2Form_Outer_Special
from objects.CSCG._2d.forms.standard._2_form.base.main import _2Form_BASE


class _2dCSCG_2Form_Outer(_2Form_BASE):
    """
    Outer standard 2-form.

    :param mesh:
    :param space:
    :param hybrid:
    :param numbering_parameters:
    :param name:
    """
    def __init__(
            self, mesh, space, hybrid=True,
            numbering_parameters='Naive',  name='outer-oriented-2-form'
    ):
        super().__init__(mesh, space, hybrid, 'outer', numbering_parameters, name)
        super().__init_2form_base__()
        self._k_ = 2
        self.standard_properties.___PRIVATE_add_tag___('2dCSCG_standard_outer_2form')
        self.standard_properties.___PRIVATE_add_tag___('2dCSCG_standard_2form')
        self._special_ = _2Form_Outer_Special(self)
        self.___kwargs___ = {
            'hybrid': hybrid,
            'numbering_parameters': numbering_parameters,
            'name': name,
        }
        self._freeze_self_()

    @property
    def special(self):
        return self._special_


if __name__ == '__main__':
    # mpiexec python _2dCSCG\form\standard\_2_form_outer.py
    # import numpy as np

    from objects.CSCG._2d.master import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.0, bounds=([0, 1], [0, 1]))([1, 1])
    # mesh = MeshGenerator('chp1',)([2,2])
    space = SpaceInvoker('polynomials')([('Lobatto', 2), ('Lobatto', 2)])
    FC = FormCaller(mesh, space)

    ES = ExactSolutionSelector(mesh)('sL:sincos1')

    f2 = FC('2-f-o', is_hybrid=False)

    # f2.TW.func.do.set_func_body_as(ES, 'potential')
    # f2.TW.current_time = 0
    # f2.TW.do.push_all_to_instant()
    # f2.discretize()
    M0 = f2.matrices.mass[0]

    print(M0.toarray())
    # save(f2, 'test_2d_f2_o')
