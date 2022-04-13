


from root.config.main import cOmm
from itertools import combinations



class OBJ_SurInternal_EdgeElement(object):
    """"""
    def __init__(self, mesh, i):
        """

        Parameters
        ----------
        mesh :
        i : int
            The edge-element #i.
        """

        if i in mesh.edge.elements:

            ee = mesh.edge.elements[i]
            positions = ee.positions

            involved_mesh_elements = list()
            for pos in positions:
                assert pos[:-2].isnumeric()
                involved_mesh_elements.append(int(pos[:-2]))

        else:
            involved_mesh_elements = None

        involved_mesh_elements = cOmm.allgather(involved_mesh_elements)

        for _ in involved_mesh_elements:
            if _ is not None:
                involved_mesh_elements = _
                break

        self._useful_ = False # this local core is useful for this connection.
        for m in involved_mesh_elements:
            if m in mesh.elements:
                self._useful_ = True
                break

        involved_mesh_elements.sort()

        COMs = combinations(involved_mesh_elements, 2)
        Ts = dict()
        for com in COMs:
            T = mesh.elements.do.find.trace_element_between_two_elements(*com)
            LT = len(T)
            if LT == 2:
                for t in T:
                    ees = mesh.trace.elements.do.find.edge_elements_surrounding_element(t)
                    if i in ees:
                        break
            elif LT == 1:
                t = T[0]
            elif LT == 0:
                t = None
            else:
                raise Exception()

            if t is not None:
                Ts[t] = com

        sequence = [involved_mesh_elements[0],]

        while len(Ts) > 0:
            start_point = sequence[-1]

            for te in Ts:
                if start_point in Ts[te]:
                    break

            sequence.append(te)

            p1, p2 = Ts[te]

            if start_point == p1:
                sequence.append(p2)
            else:
                assert start_point == p2
                sequence.append(p1)

            del Ts[te]

        assert sequence[0] == sequence[-1]

        if self._useful_:
            self._sequence_ = sequence


    @property
    def sequence(self):
        """list: A list of 9 entries represent the connection of surrounding mesh-elements and
        trace-elements.

        So `sequence[0]` must be equal to `sequence[-1]` which is a mesh-element.

        For example, for an internal edge element i whose surrounding objects are

             |  |  |
             |  |  |
        0    |  2  |     2
             |  |  |
        ------     ---------
        --4---  o  -----12--
        ------     ---------
             |  |  |
        8    |  30 |    10
             |  |  |
             |  |  |

        The sequence will be `[0, 2, 2, 12, 10, 30, 8, 4, 0]`

         m   t   m   t   m   t   m   t   m
         |   |   |   |   |   |   |   |   |
         v   v   v   v   v   v   v   v   v
        [0,  2,  2, 12, 10, 30,  8,  4,  0]

        And we only return this sequence in the cores that involve at least one of the surrounding
        mesh-element.

        """
        if self._useful_:
            return self._sequence_
        else:
            return None

