# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('../')
from screws.frozen import FrozenOnly
from _2dCSCG.forms.standard._0_form.base import _0Form_BASE

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
        self._k_ = 0
        self.standard_properties.___PRIVATE_add_tag___('2dCSCG_standard_inner_0form')
        self.standard_properties.___PRIVATE_add_tag___('2dCSCG_standard_0form')
        self._special_ = _0Form_Inner_Special(self)
        self.RESET_cache()
        self._freeze_self_()

    @property
    def special(self):
        return self._special_






class _0Form_Inner_Special(FrozenOnly):
    def __init__(self, _0sf):
        self._sf_ = _0sf
        self._freeze_self_()





if __name__ == '__main__':
    # mpiexec python _2dCSCG\form\standard\_0_form_inner.py

    from _2dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.3)([50,45])
    # mesh = MeshGenerator('chp1',)([2,2])
    space = SpaceInvoker('polynomials')([('Lobatto',3), ('Lobatto',4)])
    FC = FormCaller(mesh, space)

    ES = ExactSolutionSelector(mesh)('sL:sincos1')

    # # mesh.domain.visualize()
    #
    # def p(t, x, y): return np.cos(2*np.pi*x) * np.cos(2*np.pi*y) + t/2
    # scalar = FC('scalar', p)

    f0 = FC('0-f-i', is_hybrid=True)
    f0.TW.func.DO.set_func_body_as(ES, 'potential')
    f0.TW.current_time = 0
    f0.TW.DO.push_all_to_instant()
    f0.discretize()
    print(f0.error.L())

    f0.visualize.matplot()