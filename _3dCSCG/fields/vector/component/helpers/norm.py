

class ___VECTOR_NORM_COMPONENT___(object):
    """Here we wrap the reconstruction of standard vector such that it works like a function. Then we
    can use it to, for example, build more vectors."""
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
            w_dot_n = w[0]*n[0] + w[1]*n[1] + w[2]*n[2]
            norm_component = (w_dot_n * n[0], w_dot_n * n[1], w_dot_n * n[2])
            return xyz, norm_component

        else:
            raise Exception(f"We cannot compute norm component for {self._vf_.ftype} _3dCSCG_VectorField.")