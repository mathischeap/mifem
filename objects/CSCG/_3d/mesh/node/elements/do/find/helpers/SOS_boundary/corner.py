


from screws.freeze.base import FrozenOnly


class _3dCSCG_CornerNodeSOS(FrozenOnly):
    """"""

    def __init__(self, mesh, i):
        """"""
        trace_elements = mesh.trace.elements
        edge_elements = mesh.edge.elements
        elements = mesh.elements


        if i not in mesh.node.elements:
            self._mesh_ = None
            return

        node = mesh.node.elements[i]
        positions = node.positions
        mesh_element_corner_name = positions[0]
        mesh_element = int(mesh_element_corner_name[:-3])
        corner_name = mesh_element_corner_name[-3:]

        self._mesh_ = (mesh_element, corner_name)


        s0, s1, s2 = corner_name

        side0, corner0 = s0, s1 + s2
        side1, corner1 = s1, s0 + s2
        side2, corner2 = s2, s0 + s1

        ___ = 'NSWEBF'
        T_MAP = trace_elements.map[mesh_element]
        E_MAP = elements.map[mesh_element]
        self._trace_ = list()
        self._mesh_boundaries_ = list()
        for side, corner in zip((side0, side1, side2), (corner0, corner1, corner2)):
            _ = ___.index(side)
            trace = T_MAP[_]
            boundary = E_MAP[_]
            assert isinstance(boundary, str), f"we must find a mesh boundary."
            self._trace_.append((trace, corner, boundary))
            self._mesh_boundaries_.append(boundary)



        ___ = ['WB', 'EB', 'WF', 'EF', 'NB', 'SB', 'NF', 'SF', 'NW', 'SW', 'NE', 'SE']
        E_MAP = edge_elements.map[mesh_element]
        self._edge_ = list()
        for end, corner_edge in zip((side0, side1, side2), (corner0, corner1, corner2)):
            edge = E_MAP[___.index(corner_edge)]
            self._edge_.append((edge, end))



        self._freeze_self_()

    @property
    def mesh(self):
        return self._mesh_

    @property
    def trace(self):
        return self._trace_

    @property
    def edge(self):
        return self._edge_

    @property
    def mesh_boundaries(self):
        return self._mesh_boundaries_