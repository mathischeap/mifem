# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/04 10:40 PM
"""

import sys
if './' not in sys.path:
    sys.path.append('./')

from root.config.main import COMM, RANK, MASTER_RANK
from importlib import import_module
from components.freeze.base import FrozenOnly
from components.miscellaneous.timer import MyTimer
from components.miscellaneous.miprint import miprint
from objects.miUsGrid.triangular.mesh.main import miUsGrid_TriangularMesh
from objects.miUsGrid.triangular.space.main import miUsGrid_TriangularFunctionSpace

from objects.miUsGrid.triangular.forms.allocator import miUsGrid_FormsAllocator

from objects.miUsGrid.triangular.exactSolution.allocator import miUsGrid_ExactSolutionAllocator

from objects.miUsGrid.triangular.fields.scalar.main import miUsGrid_Triangular_Scalar
from objects.miUsGrid.triangular.fields.vector.main import miUsGrid_Triangular_Vector
from objects.miUsGrid.triangular.mesh.instances.samples.allocator import miUsGrid_TriangularMeshAllocator
from objects.miUsGrid.triangular.mesh.instances.realtime.allocator import miUsGrid_RealTime_TriangularMeshAllocator


class Call(FrozenOnly):
    """"""

    def __init__(self, mesh_source_or_ID, p,
                 boundaries=None, show_info=False, mesh_name='NoNameMesh', **kwargs):
        """

        Parameters
        ----------
        mesh_source_or_ID
        p
        boundaries
        show_info
        mesh_name
        kwargs :
            For generating real time mesh.
        """

        if mesh_source_or_ID in miUsGrid_TriangularMeshAllocator.___mesh_path___():
            mesh_name = mesh_source_or_ID
            # do not change the sequence of below two lines.
            boundaries = miUsGrid_TriangularMeshAllocator.___mesh_boundaries___()[mesh_source_or_ID]
            mesh_source_or_ID = miUsGrid_TriangularMeshAllocator.___mesh_path___()[mesh_source_or_ID]

        elif miUsGrid_RealTime_TriangularMeshAllocator.check_mesh(mesh_source_or_ID):
            mesh_name = mesh_source_or_ID
            mesh_source_or_ID, boundaries = miUsGrid_RealTime_TriangularMeshAllocator.make_mesh(
                mesh_source_or_ID, **kwargs)

        else:
            pass

        COMM.barrier()
        if show_info:
            miprint(f"---[miUsTriangle]-{MyTimer.current_time()} ...... ")
            miprint(f"   Mesh name: {mesh_name}")
        self._mesh_ = miUsGrid_TriangularMesh(mesh_source_or_ID, boundaries, name=mesh_name)
        if show_info:
            miprint(f"   Total elements: {self._mesh_.elements.num.global_cells}")
        self._space_ = miUsGrid_TriangularFunctionSpace(p)
        if show_info:
            miprint(f"   Polynomial degree: {self._space_.p}")
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
            return getattr(import_module(form_classpath), form_classname)(
                self.mesh, self.space, *args, **kwargs)

        elif ID in miUsGrid_ExactSolutionAllocator.___es_name___():
            es_classname = miUsGrid_ExactSolutionAllocator.___es_name___()[ID]
            es_classpath = miUsGrid_ExactSolutionAllocator.___es_path___()[ID]
            return getattr(import_module(es_classpath), es_classname)(
                self.mesh, *args, **kwargs)

        else:
            raise NotImplementedError(f"form id={ID} not implemented")

    @classmethod
    def listing(cls, printing=True, returning=False):
        """For an allocator class, this list all the possibilities ONLY in the master core."""
        if RANK != MASTER_RANK:
            return

        included_allocators = [
            miUsGrid_FormsAllocator,
            miUsGrid_ExactSolutionAllocator,
        ]
        listing = '\n-------- miUsFields:\n' \
                  '>>> scalar ~ miUsGrid_Triangular_Scalar\n\n' \
                  '>>> vector ~ miUsGrid_Triangular_Vector\n\n'
        for alc in included_allocators:
            listing += '\n-------- ' + alc.__name__ + ':\n'
            listing += alc.listing(printing=False, returning=True)

        if printing:
            print(listing)
        else:
            pass
        if returning:
            return listing
        else:
            pass


if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/master.py
    fc = Call('rand0', 2, show_info=True)

    fc.listing()
