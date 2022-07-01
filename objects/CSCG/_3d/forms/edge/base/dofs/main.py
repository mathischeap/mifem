# -*- coding: utf-8 -*-
import sys
if './' not in sys.path: sys.path.append('./')
from screws.freeze.main import FrozenOnly


from objects.CSCG._3d.forms.edge.base.dofs.dof.main import _3dCSCG_Edge_forms_DOF
from objects.CSCG._3d.forms.edge.base.dofs.do.main import _3dCSCG_EdgeDofs_Do

class _3dCSCG_Edge_forms_DOFs(FrozenOnly):
    """The class of all dofs of a trace form."""
    def __init__(self, ef):
        self._ef_ = ef
        self._GM_ = ef.numbering.gathering
        self._local_range_ = self._GM_.local_range
        self._do_ = None
        self._freeze_self_()

    def __iter__(self):
        for i in range(self._GM_.GLOBAL_num_dofs):
            yield i

    def __contains__(self, i):
        """If dof #i is contained in this core?"""
        if i % 1 != 0:
            return False
        else:
            if 0 <= i < self._GM_.GLOBAL_num_dofs:
                return True
            else:
                return False

    def __getitem__(self, i):
        """Return a dof object for this particular local dof #i."""
        assert i in self, f"dof#{i} is out of range."
        return _3dCSCG_Edge_forms_DOF(self, i)


    def ___PRIVATE_FIND_local_mesh_elements_and_local_indices_of_dof___(self, i):
        """Find the local mesh-element(s) and the index(es) of its local numbering of a dof i."""


        if self._local_range_ == tuple(): # I have no local mesh elements
            return list(), list()

        if i in range(*self._local_range_):
            ELEMENTS = list()
            INDICES = list()
            # note that this does not make sure `i` is contained by the core since dofs may not fully cover the range.
            for e in self._GM_: # go through all local elements
                v = self._GM_[e]
                if i in v:
                    ELEMENTS.append(e)
                    INDICES.append(v.index(i))
                else:
                    pass
            return ELEMENTS, INDICES
        else:
            return list(), list()


    @property
    def do(self):
        if self._do_ is None:
            self._do_ = _3dCSCG_EdgeDofs_Do(self)
        return self._do_