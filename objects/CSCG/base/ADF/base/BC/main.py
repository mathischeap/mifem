



from screws.freeze.base import FrozenOnly

from objects.CSCG.base.ADF.base.BC.partial_cochain.main import PartialCochain
from objects.CSCG.base.ADF.base.BC.partial_cochain.partial_dofs.main import PartialDofs



class CSCG_ADForm_BC(FrozenOnly):
    """"""
    def __init__(self, adf):
        self._adf_ = adf
        self.___PRIVATE_reset_cache___()
        self._freeze_self_()


    def ___PRIVATE_reset_cache___(self):
        self._valid_boundaries_ = None
        self._partial_cochain_ = None
        self._partial_dofs_ = None


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
            assert bn in self._adf_.mesh.boundaries.names, \
                f"boundary_names[{i}]: {bn} is not in mesh.boundaries.names: " \
                f"{self._adf_.mesh.boundaries.names}"

        self._valid_boundaries_ = valid_boundaries
        self._partial_cochain_ = None
        self._partial_dofs_ = None



    @property
    def partial_cochain(self):
        """We will interpret the BC as a PartialCochain instance which then can
        further be interpreted as data structures that can be used by,
        for example, EWC sparse matrices.
        """
        if self._partial_cochain_ is None:
            pc = PartialCochain(self._adf_)
            pc.include.boundaries(self.valid_boundaries)
            self._partial_cochain_ = pc
        return self._partial_cochain_


    @property
    def partial_dofs(self):
        """We will interpret the BC as a PartialCochain instance which then can
        further be interpreted as data structures that can be used by,
        for example, EWC sparse matrices.
        """
        if self._partial_dofs_ is None:
            pd = PartialDofs(self._adf_)
            pd.include.boundaries(self.valid_boundaries)
            self._partial_dofs_ = pd
        return self._partial_dofs_