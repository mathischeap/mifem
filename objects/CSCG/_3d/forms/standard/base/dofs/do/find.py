# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly



class _3dCSCG_SF_dofs_FIND(FrozenOnly):
    """"""
    def __init__(self, dofs):
        self._dofs_ = dofs
        self._freeze_self_()

    def dof_at_corner_of_region(self, region_name, which_corner):
        """We find the dof as a _3dCSCG_SF_DOF instance on the `which_corner` of region `region_name`.
        """

        #------ we first find the corner mesh element ------------------------------------------

        # TODO: to be continued.


    def dof_at_corner_of_mesh_element(self, i, corner_name):
        """We find the global numbering of the dof at a corner of mesh-element #i

        IMPORTANT: only in the core has mesh element i, return None in all other cores.

        Parameters
        ----------
        i : int
        corner_name : str


        Returns
        -------

        """
        assert corner_name in ['NWB', 'SWB', 'NEB', 'SEB', 'NWF', 'SWF', 'NEF', 'SEF'], \
            f"corner_name = {corner_name} is invalid."

        sf = self._dofs_._sf_

        assert sf.k in (0, 3), f"1- or 2-form has no dof at mesh element corner."

        mesh = sf.mesh

        GM = sf.numbering.gathering # make sure this is run in all cores.

        if i not in mesh.elements: return None

        GV = GM[i]

        local_numbering = sf.numbering.local[0]

        if corner_name == 'NWB':
            i, j, k = [0, 0, 0]
        elif corner_name == 'SWB':
            i, j, k = [-1, 0, 0]
        elif corner_name == 'NEB':
            i, j, k = [0, -1, 0]
        elif corner_name == 'SEB':
            i, j, k = [-1, -1, 0]
        elif corner_name == 'NWF':
            i, j, k = [0, 0, -1]
        elif corner_name == 'SWF':
            i, j, k = [-1, 0, -1]
        elif corner_name == 'NEF':
            i, j, k = [0, -1, -1]
        elif corner_name == 'SEF':
            i, j, k = [-1, -1, -1]
        else:
            raise Exception()

        return GV.full_vector[local_numbering[i, j, k]]