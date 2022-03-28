



from screws.freeze.main import FrozenOnly




class _PartialDofs_Interpretation_Globe_(FrozenOnly):
    """The class of the local dof interpretation."""
    def __init__(self, interpretation):
        self._interpretation_ = interpretation
        pd = interpretation._pd_
        self._mesh_ = pd._mesh_
        self._form_ = pd._form_
        self._dofs_ = pd._dofs_
        self._freeze_self_()


    def __iter__(self):
        """Go through all involved local element numbers."""
        for e in self._dofs_:
            yield e

    def __getitem__(self, e):
        """interpret the indicators for the involved element #e."""
        raise NotImplementedError()

