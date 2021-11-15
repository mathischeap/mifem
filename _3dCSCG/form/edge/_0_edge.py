# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('./')
from abc import ABC



from _3dCSCG.form.edge.main import _3dCSCG_Edge



class _0Edge(_3dCSCG_Edge, ABC):
    """
    Edge 0-form.

    :param mesh:
    :param space:
    :param orientation:
    :param numbering_parameters:
    :param name:
    """
    def __init__(self, mesh, space, orientation='outer',
        numbering_parameters='Naive', name='outer-oriented-0-edge-form'):
        super().__init__(mesh, space, orientation, numbering_parameters, name)
        self._k_ = 0
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_edge_0form')
        self.RESET_cache()
        self._freeze_self_()


    def RESET_cache(self):
        super().RESET_cache()





















if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\form\edge\_0_edge.py
    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.25)([5,6,7])
    space = SpaceInvoker('polynomials')([('Lobatto',2), ('Lobatto',1), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    e0 = FC('0-e')

    print(e0.standard_properties.tags)