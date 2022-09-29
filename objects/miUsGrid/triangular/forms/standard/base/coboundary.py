# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/22 6:57 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from importlib import import_module
from screws.freeze.base import FrozenOnly
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.main import EWC_SparseMatrix



class miUs_Triangular_SF_Coboundary(FrozenOnly):
    """"""

    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._incidenceMatrix_ = None
        self._freeze_self_()

    @property
    def incidence_matrix(self):
        assert self._sf_.k < 2, "volume form has no incidence matrix."
        if self._incidenceMatrix_ is None:
            formName = self._sf_.__class__.__name__
            _incidenceMatrix_ = getattr(self._sf_.space.incidence_matrix, formName)
            self._incidenceMatrix_ = \
                EWC_SparseMatrix(self._sf_.mesh.elements, _incidenceMatrix_, 'constant')
        return self._incidenceMatrix_



    def ___PRIVATE_next_class___(self):
        assert self._sf_.k < 2, "volume form has no next prime space."
        k = self._sf_.k
        base_path = '.'.join(str(self.__class__).split(' ')[1][1:-2].split('.')[:-3]) + '.'
        nextPath = base_path + f'_{k+1}.' + self._sf_.orientation + '.main'
        nextName = f'miUsTriangular_S{k+1}F_'
        if self._sf_.orientation == 'inner':
            nextName += 'Inner'
        elif self._sf_.orientation == 'outer':
            nextName += 'Outer'
        else:
            raise Exception()
        return getattr(import_module(nextPath), nextName)

    @property
    def cochain_local(self):
        """Return the local cochain (not the form) of its coboundary.

        :return:
        """
        selfCochain = self._sf_.cochain.local
        nextCochain = dict()
        incidence_matrix = self.incidence_matrix
        for i in self._sf_.mesh.elements:
            nextCochain[i] = incidence_matrix[i] @ selfCochain[i]
        return nextCochain

    def __call__(self):
        """
        When we call the coboundary object, we do the ``coboundary`` process; let ``self`` be a ``k``-form,
        it returns a ``(k+1)``-form.

        :return: A new standard ``(k+1)``-form.
        :raise AssertionError: If ``self`` has no cochain.
        """
        assert self._sf_.cochain.local is not None, "I need a cochain to perform coboundary."
        nextFmClass = self.___PRIVATE_next_class___()
        nextFmInstance = nextFmClass(
            self._sf_.mesh, self._sf_.space,
            name = 'd(' + self._sf_.standard_properties.name + ')'
        )
        nextFmInstance.cochain.local = self.cochain_local
        return nextFmInstance


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
