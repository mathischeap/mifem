
import sys
if './' not in sys.path: sys.path.append('./')
from SCREWS.frozen import FrozenOnly
from root.config import rAnk
from _3dCSCG.form.trace.dofs.basis_function import _3dCSCG_TF_DOF_BF



class _3dCSCG_Trace_forms_DOFs(FrozenOnly):
    """The class of all dofs of a trace form."""
    def __init__(self, tf):
        self._tf_ = tf
        self._GM_ = tf.numbering.gathering
        self._local_range_ = self._GM_.local_range
        self._visualize_ = None
        self._freeze_self_()


    def __iter__(self):
        assert self._tf_.IS_hybrid, f"trace form must be hybrid."
        for me in self._GM_: # go through all local mesh elements
            for i in self._GM_[me]: # go through all local dofs
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
        return _3dCSCG_Trace_forms_DOF(self, i)


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




class _3dCSCG_Trace_forms_DOF(FrozenOnly):
    """A dof of a trace form."""
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
        self._tf_ = dofs._tf_
        self._bf_ = None
        self._freeze_self_()

    @property
    def basis_function(self):
        """The local basis function(s) of this dof."""
        if self._bf_ is None:
            self._bf_ = _3dCSCG_TF_DOF_BF(self)
        return self._bf_







if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\form\trace\dofs\main.py
    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.3)(None)
    print(mesh.elements.layout)
    space = SpaceInvoker('polynomials')([('Lobatto',1), ('Lobatto',1), ('Lobatto',1)])
    FC = FormCaller(mesh, space)


    t0 = FC('0-t')

    dofs = t0.dofs
    if 0 in dofs:
        DI = dofs[0]
        bf = DI.basis_function