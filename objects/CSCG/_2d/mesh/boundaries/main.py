"""In `mesh.domain`, we distinguish domain boundary and periodic boundary.

So even two domain boundaries are periodic boundaries (physically, they are not boundaries),
we still have them named and in `mesh.domain.regions.map` shown.

However, in mesh.boundaries, the periodic boundaries are not considered. Because they have no
differences to element boundary.

Therefore, for a periodic domain, mesh.boundaries will have no valid boundary.

This is very important. The reason we have this is because of the logic we used to code the mesh. We first
generate the mesh.elements.map through `regions.map`, then we adjust the elements.map through studying the
periodic setting. This I know is not very good. But the thing is when I first code it, I did not consider
periodic boundaries. So ...

"""
import sys
if './' not in sys.path: sys.path.append('./')
from root.config.main import *
from screws.freeze.main import FrozenOnly

from objects.CSCG._2d.mesh.boundaries.visualize.main import _2dCSCG_Mesh_Boundaries_Visualize
from objects.CSCG._2d.mesh.boundaries.boundary.main import _2dCSCG_Mesh_Boundary



class _2dCSCG_Mesh_Boundaries(FrozenOnly):
    """"""
    def __init__(self, mesh):
        assert mesh.__class__.__name__ == '_2dCSCG_Mesh'
        self._mesh_ = mesh
        self._visualize_ = _2dCSCG_Mesh_Boundaries_Visualize(self)
        self._involved_elements_ = None
        self._range_of_elements_ = None
        self._boundaries_dict_ = dict()
        self._names_ = None
        self._RANGE_element_edges_ = None
        self._RANGE_trace_elements_ = None
        self._freeze_self_()

    def RESET_cache(self):
        self._involved_elements_ = None
        self._range_of_elements_ = None
        self._boundaries_dict_ = dict()
        self._names_ = None
        self._RANGE_element_edges_ = None
        self._RANGE_trace_elements_ = None

    def ___PRIVATE_parse_boundaries___(self):
        """We study the elements.map and trace.elements.map to get information we need."""
        names = list()
        Res = dict()
        Rte = dict()

        side_names = 'UDLR'
        elements = self._mesh_.elements
        t_elements = self._mesh_.trace.elements
        bns = self._mesh_.domain.boundaries.names
        for i in elements.map:
            assert i in t_elements.map, "A trivial check."
            for j in range(4):
                target = elements.map[i][j]

                if isinstance(target, str):

                    if target in bns:
                        if target not in names: names.append(target)
                        if target not in Res: Res[target] = list()
                        if target not in Rte: Rte[target] = list()
                        Res[target].append(str(i)+side_names[j])
                        Rte[target].append(t_elements.map[i][j])

                    else:
                        raise Exception("In structured mesh, str entry in elements.map "
                                        "should only be boundary name")
                else:
                    pass
        names = set(names)

        # Now, gather names ...
        names = COMM.gather(names, root=MASTER_RANK)
        if RANK == MASTER_RANK:
            NAMES = set()
            for _ in names:
                NAMES.update(_)
            names = list(NAMES)
        names = COMM.bcast(names, root=MASTER_RANK)
        # ...

        for na in names:
            if na not in Res: Res[na] = list()
            if na not in Rte: Rte[na] = list()
        self._names_ = names
        self._RANGE_element_edges_ = Res
        self._RANGE_trace_elements_ = Rte

    @property
    def names(self):
        """All boundary names (all global boundaries) we have in this mesh, periodic boundaries
        are excluded.

        Therefore, even there are some boundaries having no business with the elements in this core,
        we still include their names here.
        """
        if self._names_ is None:
            self.___PRIVATE_parse_boundaries___()
        return self._names_

    @property
    def involved_elements(self):
        """all local elements without indicating where it is."""

        if self._involved_elements_ is None:

            IE = list()

            for bn in self.range_of_elements:
                IE.extend(self.range_of_elements[bn])

            self._involved_elements_ = IE

        return self._involved_elements_

    @property
    def range_of_elements(self):
        """(dict) Return a dict that contains the local elements on each boundary."""
        if self._range_of_elements_ is None:

            RoE = dict()

            REE = self.range_of_element_edges
            for bn in self.names:
                ES = REE[bn]

                LIST = list()

                for es in ES:
                    LIST.append(int(es[:-1]))

                RoE[bn] = LIST

            self._range_of_elements_ = RoE

        return self._range_of_elements_

    @property
    def range_of_element_edges(self):
        """(dict) Return a dict that contains the local element sides on each boundary."""
        if self._RANGE_element_edges_ is None:
            self.___PRIVATE_parse_boundaries___()
        return self._RANGE_element_edges_

    @property
    def range_of_trace_elements(self):
        """(dict) Return a dict that contains the local trace elements on each boundary."""
        if self._RANGE_trace_elements_ is None:
            self.___PRIVATE_parse_boundaries___()
        return self._RANGE_trace_elements_

    @property
    def visualize(self):
        return self._visualize_

    def __getitem__(self, bn):
        if bn not in self._boundaries_dict_:
            assert bn in self.names
            self._boundaries_dict_[bn] = _2dCSCG_Mesh_Boundary(self, bn)
        return self._boundaries_dict_[bn]

    def __iter__(self):
        for name in self.names:
            yield name

    def __contains__(self, item):
        return item in self.names

    def __len__(self):
        return len(self.names)







if __name__ == "__main__":
    # mpiexec python _2dCSCG\mesh\boundaries.py
    from objects.CSCG._2d.master import MeshGenerator
    mesh = MeshGenerator('crazy', bounds=((0,3),(0,3)))([2,3])
    mesh.boundaries.___PRIVATE_parse_boundaries___()
    print(RANK, mesh.boundaries.range_of_element_edges)
    print(RANK, mesh.boundaries.range_of_trace_elements)