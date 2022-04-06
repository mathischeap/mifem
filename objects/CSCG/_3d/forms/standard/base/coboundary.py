



from objects.CSCG.base.forms.standard.coboundary import CSCG_Standard_Form_Coboundary_BASE
from importlib import import_module
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.main import EWC_SparseMatrix






class _3dCSCG_Standard_Form_Coboundary(CSCG_Standard_Form_Coboundary_BASE):
    def __init__(self, sf):
        super().__init__(sf)

    @property
    def incidence_matrix(self):
        """(scipy.sparse.csc_matrix) Return ths incidence matrix of the standard form."""
        assert self._sf_.k < 3, "volume form has no incidence matrix."
        if self._incidenceMatrix_ is None:
            formName = self._sf_.__class__.__name__
            _incidenceMatrix_ = getattr(self._sf_.space.incidence_matrix, formName)
            self._incidenceMatrix_ = \
                EWC_SparseMatrix(self._sf_.mesh.elements, _incidenceMatrix_, 'constant')

        return self._incidenceMatrix_




    def ___PRIVATE_next_class___(self):
        assert self._sf_.k < 3, "volume form has no next prime space."
        k = self._sf_.k
        base_path = '.'.join(str(self.__class__).split(' ')[1][1:-2].split('.')[:-3]) + '.'
        nextPath = base_path + f'_{k+1}s.main'
        nextName = f'_3dCSCG_{k+1}Form'
        return getattr(import_module(nextPath), nextName)