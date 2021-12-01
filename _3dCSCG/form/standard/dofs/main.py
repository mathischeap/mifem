

import sys
if './' not in sys.path: sys.path.append('./')
from SCREWS.frozen import FrozenOnly
from root.config import *

from _3dCSCG.form.standard.dofs.visualize.main import _3dCSCG_SF_DOFs_VISUALIZE, _3dCSCG_SF_DOF_VISUALIZE
from _3dCSCG.form.standard.dofs.basis_function import _3dCSCG_SF_DOF_BF

class _3dCSCG_Standard_forms_DOFs(FrozenOnly):
    """The class of all dofs of a standard form."""
    def __init__(self, sf):
        self._sf_ = sf
        self._GM_ = sf.numbering.gathering
        self._local_range_ = self._GM_.local_range
        self._visualize_ = None
        self._freeze_self_()



    def __iter__(self):
        if self._sf_.IS_hybrid:
            for me in self._GM_: # go through all local mesh elements
                for i in self._GM_[me]: # go through all local dofs
                    yield i
        else:
            SET = set()
            for me in self._GM_: # go through all local mesh elements
                v = self._GM_[me]
                SET.update(v.full_vector)
            for i in SET:
                yield i

    def __contains__(self, i):
        """If dof #i is contained in this core?"""
        if self._local_range_ == tuple(): # I have no local mesh elements
            return False

        if i in range(*self._local_range_):
            # note that this does not make sure i in contained by the core since dofs may not fully cover the range.
            for e in self._GM_: # go through all local elements
                v = self._GM_[e]
                if i in v:
                    return True # we return here!
                else:
                    pass
            return False # We have checked all gathering vectors, and find nothing.
        else:
            return False

    def __getitem__(self, i):
        """Return a dof object for this particular local dof #i."""
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



class _3dCSCG_Standard_forms_DOF(FrozenOnly):
    """A dof of a standard form."""
    def __init__(self, dofs, i):
        """"""
        # we first check if dof #i is a local dof, if not, raise Error.
        ELEMENTS, INDICES = dofs.___PRIVATE_FIND_local_mesh_elements_and_local_indices_of_dof___(i)
        assert len(ELEMENTS) > 0, f"dof #{i} is not a local dof in RANK {rAnk}."
        self._local_positions_ = list()
        for E, I in zip(ELEMENTS, INDICES):
            self._local_positions_.append((E, I))
        self._i_ = i # I am the #i dof.
        self._dofs_ = dofs
        self._sf_ = dofs._sf_
        self._visualize_ = None
        self._bf_ = None
        self._freeze_self_()


    @property
    def positions(self):
        """Return a list of tuples which represent the "LOCAL" positions. For example,
            local_positions = [[(3, 4), (4, 3), (5, 2), (6, 1), (7, 0)]]
        Then we know, this dof is at GM[3][4], GM[4][3], GM[5][2], GM[6][1], GM[7][0].

        For each item of this list, for example, (5,2) means this dof is at mesh element #5, and the local
        numbering of this dof is 2, so the third local numbering. We then can identify where is it is
        according to the degree of the space and the type (k) of the form.

        """
        return self._local_positions_

    @property
    def GLOBAL_positions(self):
        """The "GLOBAL" positions of this dof. So if it is shared by multiple cores, we return all its
        positions. The positions are indicated in the same way as local positions, see `positions`."""
        raise NotImplementedError()


    @property
    def visualize(self):
        """To visualize this dof."""
        if self._visualize_ is None:
            self._visualize_ = _3dCSCG_SF_DOF_VISUALIZE(self)
        return self._visualize_

    @property
    def basis_function(self):
        """The local basis function(s) of this dof."""
        if self._bf_ is None:
            self._bf_ = _3dCSCG_SF_DOF_BF(self)
        return self._bf_





if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\form\standard\dofs\main.py
    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.3)([3,3,3])
    space = SpaceInvoker('polynomials')([('Lobatto',1), ('Lobatto',1), ('Lobatto',1)])
    FC = FormCaller(mesh, space)


    f0 = FC('0-f', is_hybrid=False)

    # print(0 in f0.dofs)
    # print(3 in range(1,3))
    # print(f0.dofs.local_range)
    dofs = f0.dofs
    if 0 in dofs:
        DI = dofs[0]
        XYZ, IN_SITE_BF = DI.basis_function.reconstruct([-1,0,1], [-1,0,1], [-1,0,1])
        print(XYZ, IN_SITE_BF)