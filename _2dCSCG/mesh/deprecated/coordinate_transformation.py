from _3dCSCG.mesh.deprecated.coordinate_transformation.structuredMeshCTBase.base import CTBase


class CoordinateTransformation(CTBase):
    """ """
    def __init__(self, mesh):
        """ """
        super().__init__(mesh)

    @property
    def Jacobian(self):
        """ Determinant of the Jacobian matrix. """
        J = self.Jacobian_matrix
        return J[0][0] * J[1][1] - J[0][1] * J[1][0]

    @property
    def inverse_Jacobian(self):
        """ Determinant of the inverse Jacobian matrix. """
        iJ = self.inverse_Jacobian_matrix
        return iJ[0][0] * iJ[1][1] - iJ[0][1] * iJ[1][0]

    @property
    def inverse_Jacobian_matrix(self):
        """ The inverse Jacobian matrix. """
        J = self.Jacobian_matrix
        reciprocal_Jacobian = 1 / (J[0][0] * J[1][1] - J[0][1] * J[1][0])
        iJ00 = + reciprocal_Jacobian * J[1][1]
        iJ01 = - reciprocal_Jacobian * J[0][1]
        iJ10 = - reciprocal_Jacobian * J[1][0]
        iJ11 = + reciprocal_Jacobian * J[0][0]
        return ((iJ00, iJ01),
                (iJ10, iJ11))

    @staticmethod
    def ___method___(J):
        """
        If we already knwo J, we can use this method to compute the rest without
        computing J again. This saves a bit time. And we do not store J to, like,
        `self._J_` to avoid additional memory occupancy if we forget to clear it.

        Parameters
        ----------
        J :
            Jacobian matrix

        Returns
        -------
        output1:
            The inverse_Jacobian_matrix.
        Jacobian :
            Jacobian.
        reciprocal_Jacobian :
            1 / Jacobian

        """
        Jacobian = (J[0][0] * J[1][1] - J[0][1] * J[1][0])
        reciprocal_Jacobian = 1 / Jacobian
        iJ00 = + reciprocal_Jacobian * J[1][1]
        iJ01 = - reciprocal_Jacobian * J[0][1]
        iJ10 = - reciprocal_Jacobian * J[1][0]
        iJ11 = + reciprocal_Jacobian * J[0][0]
        return ((iJ00, iJ01),
                (iJ10, iJ11)), Jacobian, reciprocal_Jacobian