from screws.freeze.inheriting.frozen_only import FrozenOnly


class OnMeshElement_Standard(FrozenOnly):
    def __init__(self, tf):
        self._tf_ = tf
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

        SELF = self._tf_

        xyz = dict()
        value = dict()

        func = SELF.___DO_evaluate_func_at_time___()

        if isinstance(i, int):
            INDICES = [i,]
        elif i is None:
            INDICES = SELF.mesh.elements.indices

        else:
            raise NotImplementedError(f"_3dCSCG_TensorField of 'standard' ftype"
                                      f" mesh-element-reconstruction currently doesn't accept i={i}.")

        for i in INDICES:
            element = SELF.mesh.elements[i]
            xyz_i = element.coordinate_transformation.mapping(xi, eta, sigma)
            v00, v01, v02 = func[0][0](*xyz_i), func[0][1](*xyz_i), func[0][2](*xyz_i)
            v10, v11, v12 = func[1][0](*xyz_i), func[1][1](*xyz_i), func[1][2](*xyz_i)
            v20, v21, v22 = func[2][0](*xyz_i), func[2][1](*xyz_i), func[2][2](*xyz_i)
            if ravel:
                xyz[i] = [I.ravel('F') for I in xyz_i]
                value[i] = ([v00.ravel('F'), v01.ravel('F'), v02.ravel('F')],
                            [v10.ravel('F'), v11.ravel('F'), v12.ravel('F')],
                            [v20.ravel('F'), v21.ravel('F'), v22.ravel('F')])
            else:
                xyz[i] = xyz_i
                value[i] = ([v00, v01, v02],
                            [v10, v11, v12],
                            [v20, v21, v22])

        return xyz, value