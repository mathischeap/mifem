
from screws.freeze.base import FrozenOnly


class CSCG_standard_form_NUM(FrozenOnly):
    """"""
    def __init__(self, f):
        self._f_ = f
        self._freeze_self_()


    @property
    def basis(self):
        """(int) Return a int which represent the number of basis function one element has."""
        return self._f_._NUM_basis_

    @property
    def basis_components(self):
        """
        (Tuple[int]) Return a tuple of integers. If it is 0- or 3-form, then it is like (x,),
        and when it is 1- or 2-form, then it is like (x,y,z) because it represents a vector. But after all,
        sum(NUM_basis_components) == NUM_basis.
        """
        return self._f_._NUM_basis_components_

    @property
    def dofs(self):
        """(int) Return number of dofs in this core."""
        return self._f_.numbering.num_of_dofs_in_this_core

    @property
    def GLOBAL_dofs(self):
        return self._f_.numbering.gathering.GLOBAL_num_dofs