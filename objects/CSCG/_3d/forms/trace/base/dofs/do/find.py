
from components.freeze.base import FrozenOnly



class _3dCSCG_TraceDofs_DoFind(FrozenOnly):
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
        tf = self._dofs_._tf_

        assert tf.k == 0 , f"1- or 2-trace-form has no dof at mesh element corner."

        mesh = tf.mesh

        GM_TEW = tf.numbering.trace_element_wise # make sure this is run in all cores.

        trace_elements = mesh.trace.elements

        if i not in mesh.elements: return None

        NSWEBF = 'NSWEBF'
        MAP = trace_elements.map[i]

        DOFS = list()

        local_numbering = tf.numbering.local

        s0, s1, s2 = corner_name

        local = local_numbering[s0][0]
        I = 0 if s1 == 'W' else -1
        J = 0 if s2 == 'B' else -1
        local = local[I, J]
        DOFS.append(GM_TEW[MAP[NSWEBF.index(s0)]][local])

        local = local_numbering[s1][0]
        I = 0 if s0 == 'N' else -1
        J = 0 if s2 == 'B' else -1
        local = local[I, J]
        DOFS.append(GM_TEW[MAP[NSWEBF.index(s1)]][local])

        local = local_numbering[s2][0]
        I = 0 if s0 == 'N' else -1
        J = 0 if s1 == 'W' else -1
        local = local[I, J]
        DOFS.append(GM_TEW[MAP[NSWEBF.index(s2)]][local])

        return DOFS


