# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('./')
from objects.CSCG._2d.forms.standard._0_form.inner.special.main import _0Form_Inner_Special
from objects.CSCG._2d.forms.standard._0_form.base.main import _0Form_BASE

class _2dCSCG_0Form_Inner(_0Form_BASE):
    """
    Standard 0-form.

    :param mesh:
    :param space:
    :param is_hybrid:
    :param numbering_parameters:
    :param name:
    """
    def __init__(self, mesh, space, is_hybrid=True,
        numbering_parameters='Naive',  name='inner-oriented-0-form'):
        super().__init__(mesh, space, is_hybrid, 'inner', numbering_parameters, name)
        super().__init_0form_base__()
        self._k_ = 0
        self.standard_properties.___PRIVATE_add_tag___('2dCSCG_standard_inner_0form')
        self.standard_properties.___PRIVATE_add_tag___('2dCSCG_standard_0form')
        self._special_ = _0Form_Inner_Special(self)
        self.RESET_cache()
        self._freeze_self_()

    @property
    def special(self):
        return self._special_








if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_2d/forms/standard/_0_form/inner/main.py
    import numpy as np

    from objects.CSCG._2d.master import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.3)([5,5], show_info=True)
    # mesh = MeshGenerator('chp1',)([2,2])
    # mesh.visualize()

    space = SpaceInvoker('polynomials')([3, 4])
    FC = FormCaller(mesh, space)
    #
    # ES = ExactSolutionSelector(mesh)('sL:sincos1')
    #
    # mesh.domain.visualize()
    def p(t, x, y): return np.cos(2*np.pi*x) * np.cos(2*np.pi*y) + t/2
    scalar = FC('scalar', p)

    f0 = FC('0-f-i', is_hybrid=True)
    f0.CF = scalar
    f0.CF.current_time = 0
    f0.discretize()

    f0.visualize.matplot.contour()