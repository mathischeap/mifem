# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('/')
from abc import ABC



from objects.CSCG._3d.forms.node.base.main import _3dCSCG_Node



class _3dCSCG_0Node(_3dCSCG_Node, ABC):
    """
    Node 0-form.

    :param mesh:
    :param space:
    :param orientation:
    :param numbering_parameters:
    :param name:
    """
    def __init__(self, mesh, space, orientation='outer',
        numbering_parameters='Naive', name='outer-oriented-0-node-form'):
        super().__init__(mesh, space, orientation, numbering_parameters, name)
        self._k_ = 0
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_node_0form')
        self.___PRIVATE_reset_cache___()
        self._freeze_self_()


    def ___PRIVATE_reset_cache___(self):
        super().___PRIVATE_reset_cache___()





















if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\form\node\_0_node.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.25)([5,6,7])
    space = SpaceInvoker('polynomials')([('Lobatto',2), ('Lobatto',1), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    n0 = FC('0-n')

    print(n0.standard_properties.tags, n0.p, n0.dqp)