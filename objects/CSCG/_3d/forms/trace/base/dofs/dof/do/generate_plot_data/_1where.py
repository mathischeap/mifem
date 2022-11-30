
from components.freeze.base import FrozenOnly
from root.config.main import np


class _3dCSCG_T1F_DOF_Where(FrozenOnly):
    """"""
    def __init__(self, dof):
        self._dof_ = dof
        self._freeze_self_()

    def __call__(self, density=10, zoom=1):
        """Only return the data in the cores having the trace elements, else, return None.

        Parameters
        ----------
        density

        Returns
        -------

        """

        trace_element_position = self._dof_.trace_element_position

        if trace_element_position is None: return

        trace_element, local_numbering, component, component_indices =  trace_element_position[:4]

        mesh = self._dof_._tf_._mesh_
        trace_element =mesh.trace.elements[trace_element]

        MAPPING = trace_element.coordinate_transformation.mapping

        nodes = self._dof_._tf_.space.nodes
        direction = trace_element.normal_direction

        if direction == 'NS':
            NODES = (nodes[1], nodes[2])
        elif direction == 'WE':
            NODES = (nodes[0], nodes[2])
        elif direction == 'BF':
            NODES = (nodes[0], nodes[1])
        else:
            raise Exception()

        i, j = component_indices
        if component == 0:
            t_node0 = NODES[0][i], NODES[0][i+1]
            t_node1 = NODES[1][j]

            dof_x = np.linspace(t_node0[0], t_node0[1], density) * zoom
            dof_y = np.ones(density) * t_node1 * zoom
        elif component == 1:
            t_node0 = NODES[0][i]
            t_node1 = NODES[1][j], NODES[1][j+1]
            dof_x = np.ones(density) * t_node0 * zoom
            dof_y = np.linspace(t_node1[0], t_node1[1], density) * zoom
        else:
            raise Exception()

        DOF_xyz = list()
        DOF_xyz.append(
            MAPPING(
                    dof_x, dof_y,
                    from_element=trace_element.CHARACTERISTIC_element,
                    side = trace_element.CHARACTERISTIC_side)
        )


        SPAN = np.linspace(-1, 1, density) * zoom
        ONES = np.ones(density) * zoom
        LINES = list()
        for n in NODES[0]:
            eps = [n * ONES, SPAN]
            xyz = MAPPING(
                    *eps,
                    from_element=trace_element.CHARACTERISTIC_element,
                    side = trace_element.CHARACTERISTIC_side)
            LINES.append(xyz)

        for n in NODES[1]:
            eps = [SPAN, n * ONES]
            xyz = MAPPING(
                  *eps,
                  from_element=trace_element.CHARACTERISTIC_element,
                  side = trace_element.CHARACTERISTIC_side)
            LINES.append(xyz)

        # make sure that if a periodic trace element is always shown twice even it is only in one core.
        if trace_element.whether.on_periodic_boundary and not trace_element.whether.shared_by_cores:
            NCP = trace_element.NON_CHARACTERISTIC_position
            element= int(NCP[:-1])
            side = NCP[-1]
            for n in NODES[0]:
                eps = [n * ONES, SPAN]
                xyz = MAPPING(
                    *eps,
                    from_element=element,
                    side=side)
                LINES.append(xyz)

            for n in NODES[1]:
                eps = [SPAN, n * ONES]
                xyz = MAPPING(
                    *eps,
                    from_element=element,
                    side=side)
                LINES.append(xyz)

            DOF_xyz.append(
                MAPPING(
                        dof_x, dof_y,
                        from_element=element,
                        side = side)
            )

        RETURN = dict()
        RETURN['Trace-Element-Frame'] = LINES
        RETURN['DOF-xyz'] = DOF_xyz

        return RETURN
