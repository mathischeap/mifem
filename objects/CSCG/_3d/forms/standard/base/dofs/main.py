# -*- coding: utf-8 -*-

import sys
if './' not in sys.path: sys.path.append('./')
from components.freeze.main import FrozenOnly

from objects.CSCG._3d.forms.standard.base.dofs.visualize.main import _3dCSCG_SF_DOFs_VISUALIZE
from objects.CSCG._3d.forms.standard.base.dofs.dof.main import _3dCSCG_Standard_forms_DOF
from objects.CSCG._3d.forms.standard.base.dofs.do.main import _3dCSCG_SF_dofs_DO

class _3dCSCG_Standard_forms_DOFs(FrozenOnly):
    """The class of all dofs of a standard form."""
    def __init__(self, sf):
        self._sf_ = sf
        self._GM_ = sf.numbering.gathering
        self._local_range_ = self._GM_.local_range
        self._visualize_ = None
        self._do_ = None
        self._freeze_self_()

    def __iter__(self):
        for i in range(self.GLOBAL_num):
            yield i

    def __contains__(self, i):
        """If dof #i is contained in this core?"""
        if i % 1 != 0:
            return False
        else:
            if 0 <= i < self.GLOBAL_num:
                return True
            else:
                return False

    def __getitem__(self, i):
        """Return a dof object for this particular local dof #i."""
        assert i in self, f"dof#{i} is out of range."
        return _3dCSCG_Standard_forms_DOF(self, i)

    def ___PRIVATE_FIND_local_mesh_elements_and_local_indices_of_dof___(self, i):
        """Find the local mesh-element(s) and the index(es) of its local numbering of a dof i."""


        if self._local_range_ == tuple(): # I have no local mesh elements
            return list(), list()

        if i in range(*self._local_range_):
            ELEMENTS = list()
            INDICES = list()
            # note that this does not make sure i in contained by the core since dofs may not fully cover the range.
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
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = _3dCSCG_SF_DOFs_VISUALIZE(self)
        return self._visualize_

    @property
    def do(self):
        if self._do_ is None:
            self._do_ = _3dCSCG_SF_dofs_DO(self)
        return self._do_

    @property
    def GLOBAL_num(self):
        return self._sf_.num.global_dofs





if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\form\standard\dofs\main.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.3)([3,3,3])
    space = SpaceInvoker('polynomials')([('Lobatto',1), ('Lobatto',1), ('Lobatto',1)])
    FC = FormCaller(mesh, space)


    f0 = FC('0-f', is_hybrid=False)

    # print(0 in f0.dofs)
    # print(3 in range(1,3))
    # print(f0.dofs.local_range)
    dofs = f0.dofs

    # print(5. in dofs)

    DI = dofs[0]

    print(DI.positions)