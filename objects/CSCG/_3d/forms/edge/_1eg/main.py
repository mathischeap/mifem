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



from objects.CSCG._3d.forms.edge.base.main import _3dCSCG_Edge
from objects.CSCG._3d.forms.edge._1eg.discretize.main import _3dCSCG_Edge1Form_Discretize
from objects.CSCG._3d.forms.edge._1eg.special import _3dCSCG_1EF_Special
from objects.CSCG._3d.forms.edge._1eg.dofs.main import _3dCSCG_E1F_Dofs

class _3dCSCG_1Edge(_3dCSCG_Edge, ABC):
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
        self._discretize_ = _3dCSCG_Edge1Form_Discretize(self)
        self._special_ = None
        self._dofs_ = None
        self._freeze_self_()



    def ___Pr_check_CF___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_ScalarField':
            assert func_body.ftype in ('standard',), \
                f"3dCSCG 1edge FUNC do not accept func _3dCSCG_ScalarField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 1edge FUNC do not accept func {func_body.__class__}")

    def ___Pr_check_BC_CF___(self, func_body):
        assert func_body.mesh.domain == self.mesh.domain
        assert func_body.ndim == self.ndim == 3

        if func_body.__class__.__name__ == '_3dCSCG_ScalarField':
            assert func_body.ftype in ('standard','boundary-wise'), \
                f"3dCSCG 1edge BC do not accept func _3dCSCG_ScalarField of ftype {func_body.ftype}."
        else:
            raise Exception(f"3dCSCG 1edge BC do not accept func {func_body.__class__}")

    @property
    def discretize(self):
        return self._discretize_


    @property
    def special(self):
        if self._special_ is None:
            self._special_ = _3dCSCG_1EF_Special(self)
        return self._special_

    @property
    def dofs(self):
        if self._dofs_ is None:
            self._dofs_ = _3dCSCG_E1F_Dofs(self)
        return self._dofs_








if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\form\edge\_1_edge.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.25)([5,6,7])
    space = SpaceInvoker('polynomials')([('Lobatto',2), ('Lobatto',1), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    e1 = FC('1-e')

    print(e1.standard_properties.tags)