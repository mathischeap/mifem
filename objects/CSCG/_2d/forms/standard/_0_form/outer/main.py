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
from objects.CSCG._2d.forms.standard._0_form.outer.special.main import _0Form_Outer_Special
from objects.CSCG._2d.forms.standard._0_form.base.main import _0Form_BASE
from objects.CSCG._2d.forms.standard._0_form.outer.boundary_integrate.main import _2dCSCG_Outer_S0Form_BI


class _2dCSCG_0Form_Outer(_0Form_BASE):
    """
    Standard outer 0-form.

    :param mesh:
    :param space:
    :param hybrid:
    :param numbering_parameters:
    :param name:
    """
    def __init__(
            self, mesh, space, hybrid=True,
            numbering_parameters='Naive',  name='outer-oriented-0-form'
    ):
        super().__init__(mesh, space, hybrid, 'outer', numbering_parameters, name)
        super().__init_0form_base__()
        self._k_ = 0
        self.standard_properties.___PRIVATE_add_tag___('2dCSCG_standard_outer_0form')
        self.standard_properties.___PRIVATE_add_tag___('2dCSCG_standard_0form')
        self._special_ = _0Form_Outer_Special(self)
        self._BI_ = _2dCSCG_Outer_S0Form_BI(self)
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
    # mpiexec -n 4 python objects/CSCG/_2d/forms/standard/_0_form/outer/main.py

    from objects.CSCG._2d.master import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.1)([5, 5])
    # mesh = MeshGenerator('chp1',)([2,2])
    space = SpaceInvoker('polynomials')([('Lobatto', 3), ('Lobatto', 3)])
    FC = FormCaller(mesh, space)
    ES = ExactSolutionSelector(mesh)('sL:sincos1')
    f0 = FC('0-f-o', hybrid=True)
    f0.CF = ES.potential
    f0.CF.current_time = 0
    f0.discretize()

    import numpy as np
    x = y = np.linspace(-1, 1, 10)
    xy, v = f0.reconstruct(x, y, element_range='Down')
    # print(xy, v)
    # f0.numbering.visualize()
