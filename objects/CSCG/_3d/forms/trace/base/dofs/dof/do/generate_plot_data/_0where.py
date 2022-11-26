
from components.freeze.base import FrozenOnly



class _3dCSCG_T0F_DOF_Where(FrozenOnly):
    """"""
    def __init__(self, dof):
        self._dof_ = dof
        self._freeze_self_()

    def __call__(self, zoom=1):
        """"""
        trace_element_position = self._dof_.trace_element_position
        tf = self._dof_._dofs_._tf_
        mesh = tf.mesh
        trace_elements = mesh.trace.elements
        x, y, z = tf.space.nodes

        if trace_element_position is not None: # this dof is in this core.

            trace_element, local_numbering, _, local_index = trace_element_position[:4]
            assert _ == 0, f"0-trace only have one component!"

            TE = trace_elements[trace_element]

            normal_direction = TE.normal_direction

            if normal_direction == 'NS':
                y *= zoom
                z *= zoom
            elif normal_direction == 'WE':
                x *= zoom
                z *= zoom
            elif normal_direction == 'BF':
                x *= zoom
                y *= zoom
            else:
                raise Exception()

            x, y, z = TE.coordinate_transformation.mapping(x, y, z, parse_3_1d_eps=True)

            x = x[local_index]
            y = y[local_index]
            z = z[local_index]

            COO = {self._dof_.i: (x, y, z)}

        else:
            COO = None

        RETURN = dict()

        RETURN['COO-xyz'] = COO

        return RETURN