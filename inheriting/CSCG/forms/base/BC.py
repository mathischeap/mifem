




from tools.CSCG.partial_cochain import PartialCochain
from screws.freeze.inheriting.frozen_only import FrozenOnly




class CSCG_Form_BC(FrozenOnly):
    def __init__(self, f):
        self._f_ = f
        self._valid_boundaries_ = None # can not put it in ___PRIVATE_reset_cache___ method
        self.___PRIVATE_reset_cache___()
        self._freeze_self_()


    @property
    def body(self):
        """
        The function body.
        """
        return self._body_

    @property
    def ftype(self):
        """
        The function type.

        :return: The function type.
        :rtype: str
        """
        return self._ftype_


    def ___PRIVATE_reset_cache___(self):
        self._body_ = None
        self._ftype_ = None
        self._partial_cochain_ = None

    @property
    def valid_boundaries(self):
        return self._valid_boundaries_

    @valid_boundaries.setter
    def valid_boundaries(self, valid_boundaries):
        """This BC is valid on ``boundary_names``."""
        if isinstance(valid_boundaries, str):
            valid_boundaries = [valid_boundaries,]

        assert isinstance(valid_boundaries, (list, tuple)), \
            f"Please put boundary names in list or tuple"

        for i, bn in enumerate(valid_boundaries):
            assert bn in self._f_.mesh.boundaries.names, \
                f"boundary_names[{i}]: {bn} is not in mesh.boundaries.names: {self._f_.mesh.boundaries.names}"

        self._valid_boundaries_ = valid_boundaries
        self._partial_cochain_ = None


    @property
    def partial_cochain(self):
        """We will interpret the BC as a PartialCochain instance which then can
        further be interpreted as data structures that can be used by,
        for example, EWC sparse matrices.
        """
        if self._partial_cochain_ is None:
            pc = PartialCochain(self._f_)
            pc.include.boundaries(self.valid_boundaries)
            self._partial_cochain_ = pc
        return self._partial_cochain_