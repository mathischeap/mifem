
"""

"""
import sys
if './' not in sys.path: sys.path.append('./')


from screws.freeze.main import FrozenOnly
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.main import EWC_SparseMatrix




class _3dCSCG_Edge_Matrices(FrozenOnly):
    """"""
    def __init__(self, ef):
        """"""
        self._ef_ = ef
        self._S_ = None
        self._freeze_self_()



    @property
    def selective(self):
        """Return the selective matrix (mesh-element -> edge element)."""
        if self._S_ is None:
            k = self._ef_.k
            formName = f'_{int(k)}Edge'
            S = getattr(self._ef_.space.selective_matrix, formName)[0]
            self._S_ = \
                EWC_SparseMatrix(self._ef_.mesh.elements, S, 'constant')
        return self._S_




if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\form\edge\matrices.py
    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.25)([5,6,7])
    space = SpaceInvoker('polynomials')([('Lobatto',2), ('Lobatto',1), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    ef = FC('1-e')

    print(ef.matrices.selective[0].toarray()[6,:])