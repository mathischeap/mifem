# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/04 10:40 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from importlib import import_module
from screws.freeze.base import FrozenOnly
from objects.miUsGrid.triangular.mesh.main import miUsGrid_TriangularMesh
from objects.miUsGrid.triangular.space.main import miUsGrid_TriangularFunctionSpace
from objects.miUsGrid.triangular.forms.allocator import miUsGrid_FormsAllocator

from objects.miUsGrid.triangular.fields.scalar.main import miUsGrid_Triangular_Scalar
from objects.miUsGrid.triangular.fields.vector.main import miUsGrid_Triangular_Vector
from objects.miUsGrid.triangular.mesh.samples.allocator import miUsGrid_TriangularMeshAllocator

class FormCaller(FrozenOnly):
    """"""

    def __init__(self, mesh_source_or_ID, p, boundaries=None):
        """"""
        if mesh_source_or_ID in miUsGrid_TriangularMeshAllocator.___mesh_path___():
            boundaries = miUsGrid_TriangularMeshAllocator.___mesh_boundaries___()[mesh_source_or_ID]
            mesh_source_or_ID = miUsGrid_TriangularMeshAllocator.___mesh_path___()[mesh_source_or_ID]
        else:
            pass
        self._mesh_ = miUsGrid_TriangularMesh(mesh_source_or_ID, boundaries)
        self._space_ = miUsGrid_TriangularFunctionSpace(p)
        self._freeze_self_()

    def __repr__(self):
        return f"miUsGrid_Triangular_Form_caller@{id(self)}=mesh-{self.mesh.name}"

    @property
    def mesh(self):
        return self._mesh_

    @property
    def p(self):
        return self.space.p

    @property
    def space(self):
        return self._space_

    def __call__(self, ID, *args, **kwargs):
        """"""
        if ID == 'scalar':
            return miUsGrid_Triangular_Scalar(self.mesh, *args, **kwargs)

        elif ID == 'vector':
            return miUsGrid_Triangular_Vector(self.mesh, *args, **kwargs)

        elif ID in miUsGrid_FormsAllocator.___form_name___():
            form_classname = miUsGrid_FormsAllocator.___form_name___()[ID]
            form_classpath = miUsGrid_FormsAllocator.___form_path___()[ID]

        else:
            raise NotImplementedError(f"form id={ID} not implemented")

        return getattr(import_module(form_classpath), form_classname)(self.mesh, self.space, *args, **kwargs)




if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/master.py
    fc = FormCaller('st2', 2)

    f0 = fc('0-f-i')
    f1 = fc('1-f-i')

    # fc.mesh.visualize()

    gm0 = f0.numbering.gathering

    # for i in gm0:
    #     print(i, gm0[i].full_vector, flush=True)

    gm1 = f1.numbering.gathering
    # for i in gm1:
    #     print(i, gm1[i].full_vector, flush=True)

    fc.mesh.visualize()

    # print(gm0.local_dofs)

    from tools.linear_algebra.gathering.regular.chain_matrix.main import Chain_Gathering_Matrix

    GM = Chain_Gathering_Matrix([gm0, gm1], chain_method='sequent')


    for i in GM:
        print(i, GM[i])