# -*- coding: utf-8 -*-


from root.config.main import COMM
from itertools import combinations





class OBJ_SurBoundary_EdgeElement(object):
    """"""
    def __init__(self, mesh, i):
        """

        Parameters
        ----------
        mesh
        i : The edge-element #i.
        """
        #----------- find involved mesh_elements, boundaries and cores -------------------------
        if i in mesh.edge.elements:

            ee = mesh.edge.elements[i]
            positions = ee.positions

            involved_mesh_elements = list()
            involved_mesh_boundaries = list()
            for pos in positions:
                if pos in mesh.boundaries.names:
                    involved_mesh_boundaries.append(pos)
                else:
                    assert pos[:-2].isnumeric()
                    involved_mesh_elements.append(int(pos[:-2]))

        else:
            involved_mesh_elements = None
            involved_mesh_boundaries = None

        involved_mesh_elements = COMM.allgather(involved_mesh_elements)
        involved_mesh_boundaries = COMM.allgather(involved_mesh_boundaries)

        for _ in involved_mesh_elements:
            if _ is not None:
                involved_mesh_elements = _
                break

        for _ in involved_mesh_boundaries:
            if _ is not None:
                involved_mesh_boundaries = _
                break

        self._useful_ = False # this local core is useful for this connection.
        for m in involved_mesh_elements:
            if m in mesh.elements:
                self._useful_ = True
                break

        involved_mesh_elements.sort()

        #--------------------- make the sequence ------------------------------------------------
        if len(involved_mesh_elements) == 1:
            # this edge element must be at the mesh boundary and at a region's corner.
            trace_elements = mesh.edge.elements.do.find.trace_elements_attached_to_element(i)

            assert len(trace_elements) == 2, f"We must find two trace-elements."

            te1, te2 = trace_elements

            mesh_element = involved_mesh_elements[0]
            if mesh_element in mesh.elements:
                bd1 = mesh.trace.elements[te1].NON_CHARACTERISTIC_position
                bd2 = mesh.trace.elements[te2].NON_CHARACTERISTIC_position
                assert bd1 in involved_mesh_boundaries and bd2 in involved_mesh_boundaries
                self._sequence_ = [bd1, te1, mesh_element, te2, bd2]
            else:
                pass

        else: # more than one mesh elements involved
            trace_elements = mesh.edge.elements.do.find.trace_elements_attached_to_element(i)

            COMs = combinations(trace_elements, 2)

            Ts = dict()
            PA = list()
            for t1_t2 in COMs:
                s_m_e = mesh.trace.elements.do.find.mesh_element_shared_by_elements(*t1_t2)

                if s_m_e is not None:
                    Ts[s_m_e] = t1_t2
                    PA.extend(t1_t2)
            # we try to find the two trace-elements on the mesh boundary.
            TE_on_MB = list()
            for _ in PA:
                COUNT = PA.count(_)
                if COUNT == 1:
                    TE_on_MB.append(_)
                else:
                    assert COUNT == 2, f"An internal trace element must appear twice."
            assert len(TE_on_MB) == 2, f"Must find two trace elements on mesh boundary."

            sequence = [TE_on_MB[0],]
            while len(Ts) > 0:
                start = sequence[-1]
                for mesh_element in Ts:
                    t1_t2 = Ts[mesh_element]
                    if start in t1_t2:
                        sequence.append(mesh_element)
                        t1, t2 = t1_t2
                        if start == t1:
                            sequence.append(t2)
                        else:
                            sequence.append(t1)
                        del Ts[mesh_element]
                        break

            assert sequence[-1] == TE_on_MB[1]

            bT1, bT2 = TE_on_MB
            if bT1 in mesh.trace.elements:
                bT1 = mesh.trace.elements[bT1]
                mb1 = bT1.NON_CHARACTERISTIC_position
            else:
                mb1 = None

            if bT2 in mesh.trace.elements:
                bT2 = mesh.trace.elements[bT2]
                mb2 = bT2.NON_CHARACTERISTIC_position
            else:
                mb2 = None

            mb1 = COMM.allgather(mb1)
            mb2 = COMM.allgather(mb2)

            for _ in mb1:
                if _ is not None:
                    mb1 = _
            for _ in mb2:
                if _ is not None:
                    mb2 = _

            sequence = [mb1,] + sequence + [mb2,]

            if self._useful_:
                self._sequence_ = sequence


    @property
    def sequence(self):
        """Like the sequence of an internal edge element, but here it starts and ends with
        a mesh boundary.
        """
        if self._useful_:
            return self._sequence_
        else:
            return None