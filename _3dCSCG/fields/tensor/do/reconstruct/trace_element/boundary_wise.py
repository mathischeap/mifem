from screws.freeze.inheriting.frozen_only import FrozenOnly


class OnTraceElement_BoundaryWise(FrozenOnly):
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