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



class _1Edge(_3dCSCG_Edge, ABC):
    """
    Edge 1-form.

    :param mesh:
    :param space:
    :param orientation:
    :param numbering_parameters:
    :param name:
    """
    def __init__(self, mesh, space, orientation='outer',
        numbering_parameters='Naive', name='outer-oriented-1-edge-form'):
        super().__init__(mesh, space, orientation, numbering_parameters, name)
        self._k_ = 1
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_edge_1form')
        self.RESET_cache()
        self._freeze_self_()


    def RESET_cache(self):
        super().RESET_cache()





















if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\form\edge\_1_edge.py
    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.25)([5,6,7])
    space = SpaceInvoker('polynomials')([('Lobatto',2), ('Lobatto',1), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    e1 = FC('1-e')

    print(e1.standard_properties.tags)