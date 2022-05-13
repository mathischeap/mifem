# -*- coding: utf-8 -*-

class ___VECTOR_T_PARA_COMPONENT___(object): # Trace parallel
    """Here we wrap the reconstruction of standard vector such that it works like a function. Then we
    can use it to, for example, build more vectors.

    For parallel trace component.
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

            wXn1 = w2 * n3 - w3 * n2
            wXn2 = w3 * n1 - w1 * n3
            wXn3 = w1 * n2 - w2 * n1

            n_X_lb_wXn_rb = (n2 * wXn3 - n3 * wXn2,
                             n3 * wXn1 - n1 * wXn3,
                             n1 * wXn2 - n2 * wXn1)

            return xyz, n_X_lb_wXn_rb

        else:
            raise Exception(f"We cannot compute T_para component for {self._vf_.ftype} _3dCSCG_VectorField.")
