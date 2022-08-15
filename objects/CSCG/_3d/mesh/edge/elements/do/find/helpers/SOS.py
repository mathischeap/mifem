# -*- coding: utf-8 -*-

from screws.freeze.base import FrozenOnly


class _3dCSCG_Edge_SOS(FrozenOnly):
    """"""
    def __init__(self, mesh, i):
        """"""
        EE = mesh.edge.elements[i]
        position = EE.positions[0]
        mesh_element = int(position[:-2])
        TMAP = mesh.trace.elements.map[mesh_element]
        ind = 'NSWEBF'.index(position[-2])
        trace_element = TMAP[ind]

        self._i_ = i
        self._r_ = [mesh_element, position[-2:]]
        self._t_ = [trace_element, position[-1]]

        self._freeze_self_()

    @property
    def i(self):
        return self._i_

    @property
    def replacing(self):
        return self._r_

    @property
    def through(self):
        return self._t_