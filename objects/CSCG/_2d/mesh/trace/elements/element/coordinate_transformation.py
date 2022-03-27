

from screws.freeze.main import FrozenOnly

class _2dCSCG_Trace_Element_CoordinateTransformation(FrozenOnly):
    def __init__(self, te):
        self._te_ = te
        self._freeze_self_()


    def mapping(self, evaluation_points):
        """
        The local mapping.

        :param evaluation_points : A tuple or list of shape (ndim-1, ...).
        """
        i = self._te_.CHARACTERISTIC_element
        element_side = self._te_.CHARACTERISTIC_edge
        ep = self._te_._elements_.___generate_full_ep___(evaluation_points, element_side)
        x, y = self._te_._mesh_.elements[i].coordinate_transformation.mapping(*ep)
        return x, y

    def Jacobian_matrix(self, evaluation_points):
        """
        The local Jacobian matrix.

        :param evaluation_points : A tuple or list of shape (ndim-1, ...).
        """
        i = self._te_.CHARACTERISTIC_element
        element_edge = self._te_.CHARACTERISTIC_edge
        ep = self._te_._elements_.___generate_full_ep___(evaluation_points, element_edge)
        J = self._te_._mesh_.elements[i].coordinate_transformation.Jacobian_matrix(*ep)
        if element_edge in 'UD':
            return J[0][1], J[1][1]
        elif element_edge in 'LR':
            return J[0][0], J[1][0]
        else:
            raise Exception()

    def metric_matrix(self, ep1d):
        """ The entries of metric_matrix is normally denoted as g_{i,j}. """
        J = self.Jacobian_matrix(ep1d)
        Gi = J[0]**2 + J[1]**2
        return Gi

    def metric(self, ep1d):
        """ g, which should be det(metric_matrix): det(G). """
        return self.metric_matrix(ep1d)



