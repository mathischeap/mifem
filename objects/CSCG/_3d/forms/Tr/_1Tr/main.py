

import sys
if './' not in sys.path: sys.path.append('./')

from objects.CSCG._3d.forms.Tr.base.main import _3dCSCG_Standard_Tr
from objects.CSCG._3d.forms.Tr._1Tr.discretize.main import _3dCSCG_1Tr_Discretize

class _3dCSCG_1Tr(_3dCSCG_Standard_Tr):
    """
    Tr 1-form.

    :param mesh:
    :param space:
    :param orientation:
    :param numbering_parameters:
    :param name:
    """
    def __init__(self, mesh, space, orientation='outer',
        numbering_parameters='Naive', name='outer-oriented-1-Tr-form'):
        super().__init__(mesh, space, orientation, numbering_parameters, name)
        self._k_ = 1
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_Tr_1form')
        self.___PRIVATE_reset_cache___()
        self._discretize_ = None
        self._freeze_self_()

    def ___PRIVATE_reset_cache___(self):
        super(_3dCSCG_1Tr, self).___PRIVATE_reset_cache___()

    @property
    def discretize(self):
        if self._discretize_ is None:
            self._discretize_ = _3dCSCG_1Tr_Discretize(self)
        return self._discretize_





if __name__ == '__main__':
    # mpiexec -n 5 python objects\CSCG\_3d\forms\Tr\_1Tr\main.py

    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, \
        FormCaller  # , ExactSolutionSelector

    mesh = MeshGenerator('crazy_periodic', c=0.0)([2, 2, 2])
    space = SpaceInvoker('polynomials')([('Lobatto', 7), ('Lobatto', 8), ('Lobatto', 9)])
    FC = FormCaller(mesh, space)
    T1 = FC('1-Tr')