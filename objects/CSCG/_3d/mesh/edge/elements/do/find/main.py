"""

"""
import sys
if './' not in sys.path: sys.path.append('./')

from screws.freeze.main import FrozenOnly
from root.config.main import cOmm, MPI

from objects.CSCG._3d.mesh.edge.elements.do.find.helpers.objects_surrounding_boundary_edge_element import \
    OBJ_SurBoundary_EdgeElement
from objects.CSCG._3d.mesh.edge.elements.do.find.helpers.objects_surrounding_internal_edge_element import \
    OBJ_SurInternal_EdgeElement



class _3dCSCG_Edge_Elements_DO_FIND(FrozenOnly):
    def __init__(self, elements):
        self._elements_ = elements
        self._mesh_ = elements._mesh_
        self._freeze_self_()


    def elements_surrounding_trace_element(self, i):
        """We try to find the four edge-elements surrounding trace element #`i`.

        Parameters
        ----------
        i : int
            The trace-element.

        Returns
        -------

        """
        return self._mesh_.trace.elements.do.find.edge_elements_surrounding_element(i)


    def trace_elements_attached_to_element(self, i):
        """Find the (at least 4, at most 5) trace elements that are attached to the edge element #`i`.

        Parameters
        ----------
        i : int

        Returns
        -------

        """
        if i in self._elements_:
            element = self._elements_[i]
            positions = element.positions

            trace_elements = list()

            for pos in positions:
                if pos[:-2].isnumeric(): # this is a mesh element position
                    mesh_element = int(pos[:-2])
                    sides = pos[-2:]

                    if mesh_element in self._elements_._mesh_.elements:

                        TMAP = self._elements_._mesh_.trace.elements.map[mesh_element]

                        IND = ['NSWEBF'.index(_) for _ in sides]

                        for ind in IND:
                            trace_elements.append(TMAP[ind])

                    else:
                        pass

        else:
            trace_elements = None

        trace_elements = cOmm.allgather(trace_elements)

        TE = set()
        for te in trace_elements:
            if te is not None:
                TE.update(te)

        assert 2 <= len(TE) <= 5, f"At least two trace-elements, at most 4, attached to a edge-element."

        return list(TE)


    def objects_surrounding(self, i):
        """We try to find the objects (mesh elements and trace elements) surrounding edge element
        #i.

        Parameters
        ----------
        i

        Returns
        -------

        """
        if i in self._elements_:
            on_mesh_boundary = self._elements_[i].IS.on_mesh_boundary
        else:
            on_mesh_boundary = False

        on_mesh_boundary = cOmm.allreduce(on_mesh_boundary, op=MPI.LOR)

        if on_mesh_boundary: # on mesh boundary
            return OBJ_SurBoundary_EdgeElement(self._mesh_, i)
        else:
            return OBJ_SurInternal_EdgeElement(self._mesh_, i)














if __name__ == '__main__':
    # mpiexec -n 4 python objects\CSCG\_3d\mesh\edge\elements\do\find\main.py

    from objects.CSCG._3d.master import MeshGenerator
    elements = [1, 1, 1]
    # mesh = MeshGenerator('crazy', c=0.0, bounds=([0,3], [0,3], [0,3]))(elements)
    mesh = MeshGenerator('bridge_arch_cracked')(elements)
    edges = mesh.edge.elements

    # edges.do.find.objects_surrounding(0)
    # S = edges.do.find.objects_surrounding(1)

    # print(S.sequence)
    for i in range(edges.GLOBAL_num):
        S = edges.do.find.objects_surrounding(i)
        edges.do.illustrate_element(i)