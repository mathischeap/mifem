


from components.freeze.base import FrozenOnly
from root.config.main import COMM



class _3dCSCG_InternalNodeSOS(FrozenOnly):
    """"""
    def __init__(self, mesh, i):
        """"""
        TRACE = list()
        EDGE = list()
        N_edge = None
        S_edge = None

        if i in mesh.node.elements:

            positions = mesh.node.elements[i].positions

            for pos in positions:

                assert pos[0].isnumeric()  # must is a mesh-element position.

                if pos[-3:] == 'SEF':

                    NWB_element = int(pos[:-3])

                    if NWB_element in mesh.elements:
                        E_MAP = mesh.edge.elements.map[NWB_element]
                        N_edge = E_MAP[3]

                        T_MAP = mesh.trace.elements.map[NWB_element]
                        trace = T_MAP[1]
                        TRACE.append((trace, 'EF'))

                        # SE, SF edge
                        edge_SE = (E_MAP[11], 'F')
                        edge_SF = (E_MAP[7], 'E')
                        EDGE.extend([edge_SE, edge_SF])

                if pos[-3:] == 'SWF':

                    element = int(pos[:-3])

                    if element in mesh.elements:

                        T_MAP = mesh.trace.elements.map[element]
                        trace = T_MAP[1]
                        TRACE.append((trace, 'WF'))

                if pos[-3:] == 'SWB':

                    element = int(pos[:-3])

                    if element in mesh.elements:

                        T_MAP = mesh.trace.elements.map[element]
                        trace = T_MAP[1]
                        TRACE.append((trace, 'WB'))

                        E_MAP = mesh.edge.elements.map[element]
                        # SW, SB edge
                        edge_SW = (E_MAP[9], 'B')
                        edge_SB = (E_MAP[5], 'W')
                        EDGE.extend([edge_SW, edge_SB])

                if pos[-3:] == 'SEB':

                    element = int(pos[:-3])

                    if element in mesh.elements:

                        T_MAP = mesh.trace.elements.map[element]
                        trace = T_MAP[1]
                        TRACE.append((trace, 'EB'))

                if pos[-3:] == 'NEF':

                    SWB_element = int(pos[:-3])

                    if SWB_element in mesh.elements:
                        E_MAP = mesh.edge.elements.map[SWB_element]
                        S_edge = E_MAP[3]

        else:
            pass

        TRACE = COMM.allgather(TRACE)
        EDGE = COMM.allgather(EDGE)
        S_edge = COMM.allgather(S_edge)
        N_edge = COMM.allgather(N_edge)

        for _ in S_edge:
            if _ is not None:
                S_edge = _
                break

        for _ in N_edge:
            if _ is not None:
                N_edge = _
                break

        ___ = list()
        for _ in TRACE:
            ___.extend(_)
        TRACE = ___

        ___ = list()
        for _ in EDGE:
            ___.extend(_)
        EDGE = ___

        self._EDGE_ = EDGE
        self._TRACE_ = TRACE
        self._S_edge_ = S_edge
        self._N_edge_ = N_edge
        self._i_ = i
        self._freeze_self_()

    @property
    def i(self):
        """We are talking about node-element #`i`."""
        return self._i_

    @property
    def trace(self):
        """The trace elements on the N-S surface. The corner name represents the
        trace-element-corner on the node-element."""
        return self._TRACE_

    @property
    def edge(self):
        """The four edge-elements on the N-S surface. The side name represents the edge-element
        side on the node-element."""
        return self._EDGE_

    @property
    def S_edge(self):
        """The edge-element on the South side of the node-element"""
        return self._S_edge_

    @property
    def N_edge(self):
        """the edge-element on the North side of the node-element."""
        return self._N_edge_