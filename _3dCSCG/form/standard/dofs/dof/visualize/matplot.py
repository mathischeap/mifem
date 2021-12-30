



from SCREWS.frozen import FrozenOnly




class _3dCSCG_SF_DOF_VISUALIZE_matplot(FrozenOnly):
    """"""
    def __init__(self, dof):
        """"""
        self._dof_ = dof
        self._mesh_ = dof._sf_.mesh
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        """Default visualizer"""