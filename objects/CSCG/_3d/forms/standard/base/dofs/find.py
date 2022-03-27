from screws.freeze.base import FrozenOnly



class _3dCSCG_SF_dofs_FIND(FrozenOnly):
    """"""
    def __init__(self, dofs):
        self._dofs_ = dofs
        self._freeze_self_()

    def dof_at_corner_of_region(self, region_name, which_corner):
        """We find the dof as a _2dCSCG_SF_DOF instance on the `which_corner` of region `region_name`.
        """

        #------ we first find the corner mesh element ------------------------------------------

        # TODO: to be continued.