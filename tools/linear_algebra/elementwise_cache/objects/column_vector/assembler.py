import numpy as np
from scipy.sparse import csc_matrix

from screws.freeze.base import FrozenOnly
from tools.linear_algebra.data_structures.global_matrix.main import GlobalVector



class EWC_ColumnVector_Assembler(FrozenOnly):
    """"""
    def __init__(self, Vec):
        self._Vec_ = Vec
        self._freeze_self_()


    def __call__(self, routine=None, **kwargs):
        """Do the assembling."""
        if routine is None:
            return self.___default_routine___()
        else:
            raise Exception(f"Assembling routine = {routine} is wrong!")

    def ___default_routine___(self):
        """"""
        assert self._Vec_.gathering_matrix is not None, "I have no gathering matrix"
        GI = self._Vec_.gathering_matrix
        DEP = int(GI.GLOBAL_num_dofs)
        ROW = list()
        DAT = list()

        if len(self._Vec_) > 0:
            for i in self._Vec_:
                Vi = self._Vec_[i]
                ROW.append(GI[i][Vi.indices])
                DAT.append(Vi.data)

            ROW = np.concatenate(ROW)
            DAT = np.concatenate(DAT)
        else:
            pass

        return GlobalVector(csc_matrix((DAT, ROW, [0, len(ROW)]), shape=(DEP, 1)))