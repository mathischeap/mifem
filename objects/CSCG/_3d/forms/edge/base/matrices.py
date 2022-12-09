# -*- coding: utf-8 -*-
"""

"""
import sys
if './' not in sys.path: sys.path.append('./')

from components.freeze.main import FrozenOnly
from tools.elementwiseCache.dataStructures.objects.sparseMatrix.main import EWC_SparseMatrix




class _3dCSCG_Edge_Matrices(FrozenOnly):
    """"""
    def __init__(self, ef):
        """"""
        self._ef_ = ef
        self._S_ = None
        self._C_ = None
        self._freeze_self_()

    @property
    def selective(self):
        """Return the selective matrix (mesh-element -> edge element)."""
        if self._S_ is None:

            k = self._ef_.k
            formName = f'_3dCSCG_{int(k)}Edge'
            S = getattr(self._ef_.space.selective_matrix, formName)[0]
            self._S_ = \
                EWC_SparseMatrix(self._ef_.mesh.elements, S, 'constant')

        return self._S_

    @property
    def complement(self):
        """Return a complementary matrix (all zero) to be modified!"""
        if self._C_ is None:

            k = self._ef_.k
            num_trace_basis = getattr(self._ef_.space.num_basis, f'_3dCSCG_{k}Trace')[0][
                self._ef_.whether.hybrid
            ]
            num_edge_basis = self._ef_.num.basis
            self._C_ = EWC_SparseMatrix(
                self._ef_.mesh.elements, (num_trace_basis, num_edge_basis)
                )

        return self._C_




if __name__ == '__main__':
    # mpiexec -n 6 python objects\CSCG\_3d\forms\edge\base\matrices.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.25)([5,6,7])
    space = SpaceInvoker('polynomials')([('Lobatto',2), ('Lobatto',1), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    ef = FC('1-e')

    C = ef.matrices.complement
    print(C.shape)