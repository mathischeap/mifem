# -*- coding: utf-8 -*-
from root.config.main import *
from components.freeze.main import FrozenOnly
from objects.CSCG._3d.forms.trace.base.dofs.dof.basis_function import _3dCSCG_TF_DOF_BF

from objects.CSCG._3d.forms.trace.base.dofs.dof.visualize.main import _3dCSCG_Trace_forms_DOF_VISUALIZE
from objects.CSCG._3d.forms.trace.base.dofs.dof.do.main import _3dCSCG_TF_dof_DO
from objects.CSCG._3d.forms.trace.base.dofs.dof.whether import _3dCSCG_TraceForm_DofWhether

class _3dCSCG_Trace_forms_DOF(FrozenOnly):
    """A dof of a trace form."""
    def __init__(self, dofs, i):
        """"""
        # we first check if dof #i is a local dof, if not, raise Error.
        ELEMENTS, INDICES = dofs.___PRIVATE_FIND_local_mesh_elements_and_local_indices_of_dof___(i)
        self._local_positions_ = list()
        for E, I in zip(ELEMENTS, INDICES):
            self._local_positions_.append((E, I))
        self._i_ = i # I am the #i dof.
        self._dofs_ = dofs
        self._tf_ = dofs._tf_
        self._bf_ = None
        self._visualize_ = None
        self._GLOBAL_positions_ = None
        self._do_ = None
        self._anchor_ = None
        self._whether_ = None
        self._freeze_self_()

    @property
    def i(self):
        """I am the ith dof (I am numbered i in the gathering matrix)."""
        return self._i_

    @property
    def anchor(self):
        """An anchor is a list of one or more (x, y, z) coordinates. For a nodal dof, of course, an anchor
        is the
        coordinates of itself. For an edge (or face or volume) dof, its anchor is the mapping of
        its middle point in the reference domain.

        So, each dof have its unique anchor. And if the dof is on a periodic boundary, it could have multiple
        anchors.

        For example, if a cube is fully periodic, the corner (nonhybrid-0-trace-form) 0-dof (nodal dof) will have
        8 anchors (in an tuple)
        """
        if self._anchor_ is None:
            mesh = self._tf_.mesh
            _ = self.trace_element_position
            trace_elements = self._tf_.mesh.trace.elements

            if _ is not None:
                te, local_numbering, component, index = self.trace_element_position[:4]
                nodes = self._tf_.space.nodes

                TE = trace_elements[te]
                k = self._tf_.k
                direction = TE.normal_direction

                xyz = list()

                i0, i1 = index

                if k == 0: # this is a dof of 0-trace-form
                    raise NotImplementedError()

                elif k == 1:
                    raise NotImplementedError()

                elif k == 2:

                    if direction == 'NS':
                        nd0, nd1 = nodes[1], nodes[2]
                    elif direction == 'WE':
                        nd0, nd1 = nodes[0], nodes[2]
                    elif direction == 'BF':
                        nd0, nd1 = nodes[0], nodes[1]
                    else:
                        raise Exception()

                    assert component == 0, f"trivial check."

                    xi = (nd0[i0] + nd0[i0+1]) / 2
                    et = (nd1[i1] + nd1[i1+1]) / 2

                    if TE.whether.on_periodic_boundary:
                        positions = TE.positions
                        for pos in positions:
                            element = int(pos[:-1])
                            if element in mesh.elements:
                                x, y, z = TE.coordinate_transformation.mapping(xi, et, from_element=element)
                                x = x[0]
                                y = y[0]
                                z = z[0]
                                xyz.append(np.array([x, y, z]))

                    else:
                        x, y, z = TE.coordinate_transformation.mapping(xi, et)
                        x = x[0]
                        y = y[0]
                        z = z[0]
                        xyz.append(np.array([x, y, z]))

                else:
                    raise Exception()

                assert xyz is not None, f"must have found at least one anchor."

            else:
                xyz = None

            xyz = COMM.gather(xyz, root=MASTER_RANK)
            if RANK == MASTER_RANK:
                anchor = list()
                for xyz_i in xyz:
                    if xyz_i is not None:
                        for ac in xyz_i:
                            if len(anchor) == 0:
                                anchor.append(ac)
                            else:
                                exist = any([np.sum(np.abs(_ - ac)) < 1e-10 for _ in anchor])
                                if exist:
                                    pass
                                else:
                                    anchor.append(ac)

            else:
                anchor = None

            self._anchor_ = COMM.bcast(anchor, root=MASTER_RANK)

        return self._anchor_

    @property
    def positions(self):
        """Return a list of tuples which represent the "LOCAL" positions. For example,
            positions = [[(3, 4), (4, 3), (5, 2), (6, 1), (7, 0)]]
        Then we know, this dof is at GM[3][4], GM[4][3], GM[5][2], GM[6][1], GM[7][0]. And mesh elements
        3, 4, 5, 6, 7 are all in this core.

        For each i of this list, for example, (5,2) means this dof is at mesh element #5, and the local
        numbering of this dof is 2, so the third local numbering. We then can identify where is it is
        according to the degree of the space and the type (k) of the form.

        If
            positions = []
        Then we know this dof has is not in this core.

        """
        return self._local_positions_

    @property
    def trace_element_position(self):
        """int: this dof is on this trace element.

        If this dof is not in this core, we return None. Otherwise, we return a tuple of four
        outputs representing the location of this dof on a trace element. For example,

        (12, 5, 1, (1,2))

        This means this dof is on trace-element #12, and the trace-wise-local-numbering is 5,
        and locally, it is representing the second component (must be 1-trace form then), and the
        local component indices is (1,2).

        """
        positions = self.positions
        if positions == list():
            return None
        else:
            position = positions[0]

            nbc = self._tf_.num.basis_onside
            num_NS = nbc['N']
            num_WE = nbc['W']
            num_BF = nbc['B']

            mesh_element, local_numbering = position

            if 0 <= local_numbering < num_NS:
                side = 'N'
            else:
                local_numbering -= num_NS
                if 0 <= local_numbering < num_NS:
                    side = 'S'
                else:
                    local_numbering -= num_NS
                    if 0 <= local_numbering < num_WE:
                        side = 'W'
                    else:
                        local_numbering -= num_WE
                        if 0 <= local_numbering < num_WE:
                            side = 'E'
                        else:
                            local_numbering -= num_WE
                            if 0 <= local_numbering < num_BF:
                                side = 'B'
                            else:
                                local_numbering -= num_BF
                                assert 0 <= local_numbering < num_BF
                                side = 'F'

            T_MAP = self._tf_.mesh.trace.elements.map[mesh_element]
            index = 'NSWEBF'.index(side)
            trace_element = T_MAP[index]

            LN = self._tf_.numbering.local[side]

            _ = -1
            local_indices = None

            for _, ln in enumerate(LN):
                local_indices = np.argwhere(ln==local_numbering)
                if local_indices.shape[0] == 1:
                    break

            assert _ != -1 and local_indices is not None, f"must find it."

            return trace_element, local_numbering, _, tuple(local_indices[0])

    @property
    def global_positions(self):
        """The "GLOBAL" positions of this dof. So if it is shared by multiple cores, we return all its
        positions. The positions are indicated in the same way as the local positions, see `positions`."""
        if self._GLOBAL_positions_ is None:
            positions = self.positions
            positions = COMM.gather(positions, root=MASTER_RANK)
            if RANK == MASTER_RANK:
                GP = list()
                for PS in positions:
                    GP.extend(PS)
            else:
                GP = None
            self._GLOBAL_positions_ = COMM.bcast(GP, root=MASTER_RANK)
        return self._GLOBAL_positions_

    @property
    def basis_function(self):
        """The local basis function(s) of this dof."""
        if self._bf_ is None:
            self._bf_ = _3dCSCG_TF_DOF_BF(self)
        return self._bf_

    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = _3dCSCG_Trace_forms_DOF_VISUALIZE(self)
        return self._visualize_

    @property
    def do(self):
        if self._do_ is None:
            self._do_ = _3dCSCG_TF_dof_DO(self)
        return self._do_

    @property
    def whether(self):
        if self._whether_ is None:
            self._whether_ = _3dCSCG_TraceForm_DofWhether(self)
        return self._whether_