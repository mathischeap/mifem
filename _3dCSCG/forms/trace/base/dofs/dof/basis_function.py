"""The class for the basis function of a trace dof (not dofs)."""





from screws.frozen import FrozenOnly






class _3dCSCG_TF_DOF_BF(FrozenOnly):
    """"""
    def __init__(self, dof):
        self._dof_ = dof
        self._tf_ = dof._tf_
        self._mesh_ = self._tf_.mesh
        self._space_ = self._tf_.space
        self._freeze_self_()