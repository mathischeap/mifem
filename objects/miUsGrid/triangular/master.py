# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/04 10:40 PM
"""

import sys

if './' not in sys.path: sys.path.append('./')

from root.config.main import cOmm
from importlib import import_module
from screws.freeze.base import FrozenOnly
from objects.miUsGrid.triangular.mesh.main import miUsGrid_TriangularMesh
from objects.miUsGrid.triangular.space.main import miUsGrid_TriangularFunctionSpace

from objects.miUsGrid.triangular.forms.allocator import miUsGrid_FormsAllocator

from objects.miUsGrid.triangular.exact_solution.allocator import miUsGrid_ExactSolutionAllocator

from objects.miUsGrid.triangular.fields.scalar.main import miUsGrid_Triangular_Scalar
from objects.miUsGrid.triangular.fields.vector.main import miUsGrid_Triangular_Vector
from objects.miUsGrid.triangular.__test__.mesh_samples.allocator import miUsGrid_TriangularMeshAllocator
from objects.miUsGrid.triangular.__test__.realtime_meshes.allocator import miUsGrid_RealTime_TriangularMeshAllocator




class Call(FrozenOnly):
    """"""

    def __init__(self, mesh_source_or_ID, p, boundaries=None, **kwargs):
        """

        Parameters
        ----------
        mesh_source_or_ID
        p
        boundaries
        kwargs :
            For generating real time mesh.
        """

        if mesh_source_or_ID in miUsGrid_TriangularMeshAllocator.___mesh_path___():
            # do not change the sequence of below two lines.
            boundaries = miUsGrid_TriangularMeshAllocator.___mesh_boundaries___()[mesh_source_or_ID]
            mesh_source_or_ID = miUsGrid_TriangularMeshAllocator.___mesh_path___()[mesh_source_or_ID]

        elif miUsGrid_RealTime_TriangularMeshAllocator.check_mesh(mesh_source_or_ID):

            mesh_source_or_ID, boundaries = miUsGrid_RealTime_TriangularMeshAllocator.make_mesh(
                mesh_source_or_ID, **kwargs)

        else:
            pass

        cOmm.barrier()
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
            return getattr(import_module(form_classpath), form_classname)(self.mesh, self.space,
                                                                          *args, **kwargs)

        elif ID in miUsGrid_ExactSolutionAllocator.___es_name___():
            es_classname = miUsGrid_ExactSolutionAllocator.___es_name___()[ID]
            es_classpath = miUsGrid_ExactSolutionAllocator.___es_path___()[ID]

            return getattr(import_module(es_classpath), es_classname)(self.mesh, *args, **kwargs)

        else:
            raise NotImplementedError(f"form id={ID} not implemented")







if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/master.py
    fc = Call('rand0', 2)

    f0 = fc('0-f-i')
    f1 = fc('1-f-i')

    # fc.mesh.visualize()

    gm0 = f0.numbering.gathering

    # for i in gm0:
    #     print(i, gm0[i].full_vector, flush=True)

    gm1 = f1.numbering.gathering
    # for i in gm1:
    #     print(i, gm1[i].full_vector, flush=True)

    # fc.mesh.visualize()

    # print(gm0.local_dofs)

    from tools.linear_algebra.gathering.regular.chain_matrix.main import Chain_Gathering_Matrix

    GM = Chain_Gathering_Matrix([gm0, gm1], chain_method='sequent')

    print(GM.GLOBAL_num_dofs)

    # for i in GM:
    #     print(i, GM[i])