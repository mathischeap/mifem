

class ___VECTOR_T_PERP_COMPONENT___(object): # Trace Perpendicular
    """Here we wrap the reconstruction of standard vector such that it works like a function. Then we
    can use it to, for example, build more vectors.

    For perpendicular trace component.
    """
    def __init__(self, vf, i):
        """

        :param vf: the vector
        :param i: for trace element #i
        """
        self._vf_ = vf
        self._i_ = i
        self._te_ = vf.mesh.trace.elements[i]

    def __call__(self, t, xi, et, sg):
        """This actually includes a reconstruction. So we will call from xi, et, sg.

        :param t: at time t
        :param xi: must be a 1d array.
        :param et: must be a 1d array.
        :param sg: must be a 1d array.
        :return:
        """
        vf = self._vf_
        i = self._i_

        if vf.ftype == 'standard':

            xyz, w = vf.reconstruct(xi, et, sg,
                                    time=t,
                                    i=i,
                                    ravel=False,
                                    where='trace-element')
            xyz = xyz[i]
            w = w[i]
            n = self._te_.coordinate_transformation.unit_normal_vector(xi, et, sg, parse_3_1d_eps=True)
            w1, w2, w3 = w
            n1, n2, n3 = n

            # #------- option 1 ----------------------------------------------------------------------
            wXn = (w2 * n3 - w3 * n2,
                   w3 * n1 - w1 * n3,
                   w1 * n2 - w2 * n1)
            #--------- option 2: TEST FOR ORTHOGONAL MESH ONLY ------------------------------------
            # side = self._te_.CHARACTERISTIC_side
            # if side in 'NS':
            #     wXn =  (w1, w3, -w2)
            # elif side in 'WE':
            #     wXn =  (-w3, w2, w1)
            # elif side in 'BF':
            #     wXn =  (w2, -w1, w3)
            # else:
            #     raise Exception

            #=========================================================================================


            return xyz, wXn

        else:
            raise Exception(f"We cannot compute T_perp component for {self._vf_.ftype} _3dCSCG_VectorField.")