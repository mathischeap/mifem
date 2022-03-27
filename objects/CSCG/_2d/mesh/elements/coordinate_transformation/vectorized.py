

from screws.freeze.base import FrozenOnly
import numpy as np



class _2dCSCG_MeshElements_CT_Vectorized(FrozenOnly):
    """We will compute the results all together and put the results in higher dimensional arrays.
    The first index always refers to the mesh-elements.
    """
    def __init__(self, elements):
        self._elements_ = elements
        self._freeze_self_()

    def mapping(self, xi, et):
        """"""
        if len(self._elements_) == 0:
            return None
        mapping = list()
        for i in self._elements_:
            mapping.append(self._elements_[i].coordinate_transformation.mapping(xi, et))
        return np.array(mapping)


    def Jacobian_matrix(self, xi, et):
        """"""
        if len(self._elements_) == 0:
            return None

        if self._elements_.IS.homogeneous_according_to_types_wrt_metric:

            i = self._elements_.indices[0]
            Jacobian_matrix = self._elements_[i].coordinate_transformation.Jacobian_matrix(xi, et)
            return Jacobian_matrix

        else:
            J00, J01, J10, J11 = list(), list(), list(), list()

            J = self._elements_.coordinate_transformation.Jacobian_matrix(xi, et)

            for i in self._elements_:
                J00.append(J[i][0][0])
                J01.append(J[i][0][1])
                J10.append(J[i][1][0])
                J11.append(J[i][1][1])

            J00 = np.array(J00)
            J01 = np.array(J01)
            J10 = np.array(J10)
            J11 = np.array(J11)

            return ([J00, J01],
                    [J10, J11])



    def Jacobian(self, xi, et):
        """"""
        if len(self._elements_) == 0:
            return None

        if self._elements_.IS.homogeneous_according_to_types_wrt_metric:

            i = self._elements_.indices[0]
            Jacobian = self._elements_[i].coordinate_transformation.Jacobian(xi, et)
            return Jacobian

        else:
            Jacobian = list()

            J = self._elements_.coordinate_transformation.Jacobian(xi, et)

            for i in self._elements_:
                Jacobian.append(J[i])
            return np.array(Jacobian)

    def inverse_Jacobian_matrix(self, xi, et):
        """"""
        if len(self._elements_) == 0:
            return None

        if self._elements_.IS.homogeneous_according_to_types_wrt_metric:

            i = self._elements_.indices[0]
            inverse_Jacobian_matrix = self._elements_[i].coordinate_transformation.inverse_Jacobian_matrix(xi, et)
            return inverse_Jacobian_matrix

        else:
            iJ00, iJ01, iJ10, iJ11 = list(), list(), list(), list()

            iJ = self._elements_.coordinate_transformation.inverse_Jacobian_matrix(xi, et)

            for i in self._elements_:
                iJ00.append(iJ[i][0][0])
                iJ01.append(iJ[i][0][1])
                iJ10.append(iJ[i][1][0])
                iJ11.append(iJ[i][1][1])

            iJ00 = np.array(iJ00)
            iJ01 = np.array(iJ01)
            iJ10 = np.array(iJ10)
            iJ11 = np.array(iJ11)

            return ([iJ00, iJ01],
                    [iJ10, iJ11])


    def inverse_Jacobian(self, xi, et):
        """"""
        if len(self._elements_) == 0:
            return None

        if self._elements_.IS.homogeneous_according_to_types_wrt_metric:

            i = self._elements_.indices[0]
            inverse_Jacobian = self._elements_[i].coordinate_transformation.inverse_Jacobian(xi, et)
            return inverse_Jacobian

        else:
            inverse_Jacobian = list()

            iJ = self._elements_.coordinate_transformation.inverse_Jacobian(xi, et)

            for i in self._elements_:
                inverse_Jacobian.append(iJ[i])
            return np.array(inverse_Jacobian)