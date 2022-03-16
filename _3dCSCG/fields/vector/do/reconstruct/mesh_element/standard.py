from screws.freeze.inheriting.frozen_only import FrozenOnly


class OnMeshElement_Standard(FrozenOnly):
    def __init__(self, vf):
        self._vf_ = vf
        self._freeze_self_()

    def __call__(self, xi, eta, sigma, ravel, i):
        """

        :param xi:
        :param eta:
        :param ravel:
        :param i:
            1) self.ftype == 'standard':
                Do the reconstruction in mesh element #i. When it is None, it means all local mesh
                elements.
        :return:
        """

        SELF = self._vf_

        xyz = dict()
        value = dict()

        func = SELF.___DO_evaluate_func_at_time___()

        if isinstance(i, int):
            INDICES = [i, ]
        elif i is None:
            INDICES = SELF.mesh.elements.indices

        else:
            raise NotImplementedError(f"_3dCSCG_VectorField of 'standard' ftype"
                                      f" mesh-element-reconstruction currently doesn't accept i={i}.")

        for i in INDICES:
            element = SELF.mesh.elements[i]
            xyz_i = element.coordinate_transformation.mapping(xi, eta, sigma)
            vx_i = func[0](*xyz_i)
            vy_i = func[1](*xyz_i)
            vz_i = func[2](*xyz_i)

            if ravel:
                xyz[i] = [I.ravel('F') for I in xyz_i]
                value[i] = [vx_i.ravel('F'), vy_i.ravel('F'), vz_i.ravel('F')]
            else:
                xyz[i] = xyz_i
                value[i] = [vx_i, vy_i, vz_i]


        return xyz, value
