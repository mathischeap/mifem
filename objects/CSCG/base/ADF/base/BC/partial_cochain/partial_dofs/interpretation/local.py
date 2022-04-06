


from screws.freeze.main import FrozenOnly


class _ADF_PartialDofs_Interpretation_Local_(FrozenOnly):
    """The class of the local dof interpretation."""
    def __init__(self, interpretation):
        self._interpretation_ = interpretation
        pd = interpretation._pd_
        self._mesh_ = pd._mesh_
        self._adf_ = pd._adf_
        self._dofs_ = pd._dofs_
        self._freeze_self_()

    def __iter__(self):
        """Go through all involved local element numbers."""
        for e in self._dofs_:
            yield e

    def __len__(self):
        """How many involved local elements?"""
        return len(self._dofs_)

    def __contains__(self, e):
        """If element #e is involved?"""
        return e in self._dofs_

    def __getitem__(self, e):
        """interpret the indicators for the involved element #e."""
        indicators = self._dofs_[e]
        DOFs = list()
        for indi in indicators:
            if indi[:2] == '1-': # type-1 indicators
                side = indi[2]
                # so the numbering property must have this method implemented.
                if  self._mesh_.ndim == 3:
                    dofs = self._adf_.prime.numbering.do.\
                        find.local_dofs_on_element_side(side)
                elif self._mesh_.ndim == 2:
                    dofs = self._adf_.prime.numbering.do.\
                        find.local_dofs_on_element_edge(side)
                else:
                    raise Exception()
                DOFs.extend(dofs)
            else:
                raise NotImplementedError(
                    f"Cannot interpret as local_dofs for indicator "
                    f"{indi} for element #{e}.")
        return DOFs
