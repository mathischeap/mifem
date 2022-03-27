
import sys
if './' not in sys.path: sys.path.append('/')
from screws.freeze.main import FrozenOnly

from objects.CSCG._3d.forms.trace.base.dofs.dof.main import _3dCSCG_Trace_forms_DOF


class _3dCSCG_Trace_forms_DOFs(FrozenOnly):
    """The class of all dofs of a trace form."""
    def __init__(self, tf):
        self._tf_ = tf
        self._GM_ = tf.numbering.gathering
        self._local_range_ = self._GM_.local_range
        self._visualize_ = None
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






if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\form\trace\dofs\main.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.3)(None)
    space = SpaceInvoker('polynomials')([('Lobatto',1), ('Lobatto',1), ('Lobatto',1)])
    FC = FormCaller(mesh, space)


    t0 = FC('0-t')

    dofs = t0.dofs
    DI = dofs[0]
    bf = DI.basis_function