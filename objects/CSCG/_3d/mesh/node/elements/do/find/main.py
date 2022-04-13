

import sys
if './' not in sys.path: sys.path.append('./')
from screws.freeze.base import FrozenOnly
from root.config.main import cOmm, MPI

from objects.CSCG._3d.mesh.node.elements.do.find.helpers.SOS_internal import _3dCSCG_InternalNodeSOS

class _3dCSCG_NodeMesh_DoFind(FrozenOnly):
    """"""
    def __init__(self, elements):
        self._elements_ = elements
        self._mesh_ = elements._mesh_
        _ = self._mesh_.edge.elements # to make sure the edge.elements.map is generated!
        self._freeze_self_()


    def hybrid_singularity_overcoming_setting(self, i):
        """

        Parameters
        ----------
        i

        Returns
        -------

        """
        if i in self._elements_ and self._elements_[i].IS.on_mesh_boundary:
            omb = True
        else:
            omb = False

        omb = cOmm.allreduce(omb, op=MPI.LOR)

        if omb:
            sos = None
        else:
            sos = _3dCSCG_InternalNodeSOS(self._mesh_, i)

        if i in self._elements_:
            return sos
        else:
            return None


    def edge_elements_attached_to_element(self, i):
        """Find the edge elements that are attached to a node element #`i`.

        For example, for a typical internal node element, there will be 6 edge elements attached
        to it.

        Parameters
        ----------
        i : int
            The node element #`i`.

        Returns
        -------

        """
        if i in self._elements_:
            element = self._elements_[i]
            positions = element.positions
            EDGE_ELEMENTS = set()
            for pos in positions:
                if pos[:-3].isnumeric():
                    mesh_element = int(pos[:-3])
                    if mesh_element in self._mesh_.elements:
                        sNS, sWE, sBF = pos[-3:]

                        edge1 = sNS + sWE
                        edge2 = sNS + sBF
                        edge3 = sWE + sBF

                        indices = [
                        ['WB', 'EB', 'WF', 'EF', 'NB', 'SB',
                         'NF', 'SF', 'NW', 'SW', 'NE', 'SE'].index(_) for _ in [edge1, edge2, edge3] ]

                        E_MAP = self._mesh_.edge.elements.map[mesh_element]

                        EDGE_ELEMENTS.update([E_MAP[_] for _ in indices])
                    else:
                        pass
                else: # a boundary location.
                    pass

        else:
            EDGE_ELEMENTS = None

        EDGE_ELEMENTS = cOmm.allgather(EDGE_ELEMENTS)

        ___ = set()

        for EE in EDGE_ELEMENTS:
            ___.update(EE)

        return list(___)


    def trace_elements_attached_to_element(self, i):
        """Find the trace elements that are attached to a node element #`i`.

        For example, for a typical internal node element, there will be 12 trace elements attached
        to it.

        Parameters
        ----------
        i : int
            The node element #`i`.

        Returns
        -------

        """
        if i in self._elements_:
            element = self._elements_[i]
            positions = element.positions
            TRACE_ELEMENTS = set()
            for pos in positions:
                if pos[:-3].isnumeric():
                    mesh_element = int(pos[:-3])
                    if mesh_element in self._mesh_.elements:
                        sides = pos[-3:]

                        T_MAP = self._mesh_.trace.elements.map[mesh_element]

                        indices = ['NSWEBF'.index(_) for _ in sides]

                        TRACE_ELEMENTS.update([T_MAP[_] for _ in indices])
                    else:
                        pass
                else: # a boundary location.
                    pass

        else:
            TRACE_ELEMENTS = None

        TRACE_ELEMENTS = cOmm.allgather(TRACE_ELEMENTS)

        ___ = set()

        for TE in TRACE_ELEMENTS:
            ___.update(TE)

        return list(___)




if __name__ == '__main__':
    # mpiexec -n 4 python objects\CSCG\_3d\mesh\node\elements\do\find\main.py

    from objects.CSCG._3d.master import MeshGenerator
    elements = [2, 2, 2]
    mesh = MeshGenerator('crazy', c=0.0, bounds=([0,3], [0,3], [0,3]))(elements)
    # mesh = MeshGenerator('bridge_arch_cracked')(elements)
    nodes = mesh.node.elements

    for i in range(nodes.GLOBAL_num):
        S = nodes.do.find.hybrid_singularity_overcoming_setting(i)
        # nodes.do.illustrate_element(i)