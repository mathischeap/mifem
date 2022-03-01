"""
In `mesh.domain`, we distinguish domain boundary and periodic boundary. So even two domain boundaries
are periodic boundaries (physically, they are not boundaries), we still have them named and in
mesh.domain.regions.map shown.

However, in mesh.boundaries, the periodic boundaries are not considered.

Therefore, for a periodic domain, mesh.boundaries will have no valid boundary. But mesh.domain.boundaries has.

This is very important. The reason we have this is because of the logic we used to code the mesh. We first
generate the mesh.elements.map through `regions.map`, then we adjust the elements.map through studying the
periodic setting. This I know is not very good. But the thing is when I first code it, I did not consider
periodic boundaries. So ...
"""

import sys
if './' not in sys.path: sys.path.append('../')

from root.config import *
from screws.frozen import FrozenOnly

from _3dCSCG.mesh.boundaries.boundary.main import _3dCSCG_Mesh_Boundary




class _3dCSCG_Mesh_Boundaries(FrozenOnly):
    """"""
    def __init__(self, mesh):
        assert mesh.__class__.__name__ == '_3dCSCG_Mesh'
        self._mesh_ = mesh
        self.___PRIVATE_reset_cache___()
        self._freeze_self_()

    def ___PRIVATE_reset_cache___(self):
        self._boundaries_dict_ = None
        self._names_ = None
        self._RANGE_element_sides_ = None
        self._RANGE_trace_elements_ = None

    def ___PRIVATE_parse_boundaries___(self):
        """We study the elements.map and trace.elements.map to get information we need."""
        names = list()
        Res = dict()
        Rte = dict()

        side_names = 'NSWEBF'
        elements = self._mesh_.elements
        t_elements = self._mesh_.trace.elements
        bns = self._mesh_.domain.boundaries.names
        for i in elements.map:
            assert i in t_elements.map, "A trivial check."
            for j in range(6):
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
        names = cOmm.gather(names, root=mAster_rank)
        if rAnk == mAster_rank:
            NAMES = set()
            for _ in names:
                NAMES.update(_)
            names = tuple(NAMES)
        names = cOmm.bcast(names, root=mAster_rank)
        # ...

        for na in names:
            if na not in Res: Res[na] = list()
            if na not in Rte: Rte[na] = list()

        # note that self._names_ != mesh.domain.boundaries.names unless the domain has no periodic boundaries.
        self._names_ = names
        self._RANGE_element_sides_ = Res
        self._RANGE_trace_elements_ = Rte

    @property
    def names(self):
        """All boundary names (all global boundaries) we have in this mesh, periodic boundaries are excluded.

        Therefore, even there are some boundaries having no business with the elements in this core,
        we still include their names here.

        Therefore, it != mesh.domain.boundaries.names unless the domain has no periodic boundaries.
        """
        if self._names_ is None:
            self.___PRIVATE_parse_boundaries___()
        return self._names_

    @property
    def range_of_element_sides(self):
        """(dict) Return a dict that contains the local element sides on each boundary."""
        if self._RANGE_element_sides_ is None:
            self.___PRIVATE_parse_boundaries___()
        return self._RANGE_element_sides_

    @property
    def range_of_trace_elements(self):
        """(dict) Return a dict that contains the local trace elements on each boundary."""
        if self._RANGE_trace_elements_ is None:
            self.___PRIVATE_parse_boundaries___()
        return self._RANGE_trace_elements_


    def __getitem__(self, bn):
        if bn not in self._boundaries_dict_:
            assert bn in self.names, f"I have no boundary named {bn}."
            self._boundaries_dict_[bn] = _3dCSCG_Mesh_Boundary(self, bn)
        return self._boundaries_dict_[bn]

    def __iter__(self):
        for name in self.names:
            yield name

    def __contains__(self, item):
        return item in self.names

    def __len__(self):
        return len(self.names)






if __name__ == "__main__":
    # mpiexec python _3dCSCG\mesh\boundaries.py
    from _3dCSCG.main import MeshGenerator
    mesh = MeshGenerator('crazy')([4,2,[1,2,4,2,1]])
    # mesh = MeshGenerator('crazy', bounds=((0,3),(0,3),(0,3)))([2,1,1])
    # mesh.boundaries.___PRIVATE_parse_boundaries___()
    print(rAnk, mesh.boundaries.range_of_element_sides)
    # print(rAnk, mesh.boundaries.names)
    # print(rAnk, mesh.boundaries.names, mesh.domain.boundaries.names)