from components.freeze.base import FrozenOnly
from root.config.main import COMM, MASTER_RANK, RANK






class _2dCSCG_SF_dofs_FIND(FrozenOnly):
    """"""
    def __init__(self, dofs):
        self._dofs_ = dofs
        self._freeze_self_()

    def dof_at_corner_of_region(self, region_name, which_corner):
        """We find the dof as a _2dCSCG_SF_DOF instance on the `which_corner` of region `region_name`.
        """

        corner_element = self._dofs_._sf_.mesh.elements.do.find.element_at_region_corner(region_name, which_corner)
        GM = self._dofs_._sf_.numbering.gathering
        if corner_element not in self._dofs_._sf_.mesh.elements:
            dof_num = None
        else:
            if which_corner[0] in 'LR':
                which_corner = which_corner[::-1]

            id0 = 0 if which_corner[0] == 'U' else -1
            id1 = 0 if which_corner[0] == 'L' else -1

            dof_local_numbering = self._dofs_._sf_.numbering.local[0][id0, id1] # local numbering
            dof_num = GM[corner_element][dof_local_numbering]

        dof_num = COMM.gather(dof_num, root=MASTER_RANK)

        if RANK == MASTER_RANK:
            dof_num = [ i for i in dof_num if i is not None ]
            assert len(dof_num) == 1, f"We must only find dof as we first will only find one mesh-element"
            dof_num = dof_num[0]

        dof_num = COMM.bcast(dof_num, root=MASTER_RANK)

        return self._dofs_[dof_num]
