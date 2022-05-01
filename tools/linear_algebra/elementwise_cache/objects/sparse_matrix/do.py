

from screws.freeze.base import FrozenOnly


class EWC_SparseMatrix_Do(FrozenOnly):
    """"""
    def __init__(self, MAT):
        """"""
        self._MAT_ = MAT
        self.___locker___ = False # By default, the locker is off!
        self.___sparsity_locker___ = False # # By default, the sparsity locker is off!
        self._freeze_self_()

    def lock_assembled_matrix(self):
        """We will no longer be able to make any customization (adjust) to the MAX. And we will save
        the assembled matrix once we assemble it.
        """
        self.___locker___ = True


    def unlock_assembled_matrix(self):
        """We will be able to make any customization (adjust) to the MAX. And we will not save
        the assembled matrix once we assemble it, and if we already save it, we just clear it.
        """
        self.___locker___ = False
        self._MAT_.assembler.___AMC___ = None # clear the assembled matrix cache


    def lock_sparsity(self):
        """We lock the sparsity of local matrices (as well as the global matrix), so the sparsity
         cannot be changed anymore although the entries can still be changed.
         """
        self.___sparsity_locker___ = True


    def unlock_sparsity(self):
        """We now can change the sparsity. And we clear the assembler cache because which probably
        will be invalid when we have changed the sparsity.
        """
        self.___sparsity_locker___ = False
        self._MAT_.assembler._cache_ = None # clear the assembler cache


    def save_memory_by_cleaning_duplicate_values(self):
        """If EWC[i] == EWC[j] but EWC[i] is not EWC[j], we clean EWC[j] and make EWC[j] point EWC[i].

        Returns
        -------

        """
        # TODO: to be continued.

