from screws.freeze.base import FrozenOnly




class _3dCSCG_EdgeDofs_DoFIND(FrozenOnly):
    """"""
    def __init__(self, dofs):
        self._dofs_ = dofs
        self._freeze_self_()

    def dof_at_corner_of_mesh_element(self, i, corner_name):
        """We find the global numbering of the (3) dofs at a corner of mesh-element #i

        IMPORTANT: Return only in the core having mesh element i, return None in all other cores.

        Parameters
        ----------
        i : int
        corner_name : str

        Returns
        -------

        """
        assert corner_name in ['NWB', 'SWB', 'NEB', 'SEB', 'NWF', 'SWF', 'NEF', 'SEF'], \
            f"corner_name = {corner_name} is invalid."
        ef = self._dofs_._ef_

        assert ef.k == 0 , f"1-edge-form has no dof at mesh element corner."
        mesh = ef.mesh

        GM_EEW = ef.numbering.edge_element_wise # make sure this is run in all cores.

        edge_elements = mesh.edge.elements

        if i not in mesh.elements: return None

        DOFS = list()
        edge_names = ['WB', 'EB', 'WF', 'EF', 'NB', 'SB', 'NF', 'SF', 'NW', 'SW', 'NE', 'SE']

        s0, s1, s2 = corner_name

        edge0, side0 = s0+s1, s2
        edge1, side1 = s0+s2, s1
        edge2, side2 = s1+s2, s0

        MAP = edge_elements.map[i]

        for edge, side in zip((edge0, edge1, edge2), (side0, side1, side2)):
            index = edge_names.index(edge)
            edge_element_number = MAP[index]
            numbering = GM_EEW[edge_element_number]
            index = 0 if side in 'NWB' else -1
            DOFS.append(numbering[index])

        return DOFS