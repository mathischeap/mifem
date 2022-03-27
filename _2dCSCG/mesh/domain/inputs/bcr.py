


from _2dCSCG.mesh.domain.inputs.base import DomainInputBase
from screws.numerical._1d_space.derivative import NumericalDerivative_fx

import numpy as np

from screws.decorators.classproperty.main import classproperty


def X(t): return 5 * t


def Xt(t): return 5 + 0 * t


def Y(t):
    return -0.3 * np.sin(np.pi * 5 * t) * (1 - np.abs(5 * t - 2.5) / 2.5)


def Yt1(t):
    return -np.pi * 5 * 0.3 * np.cos(np.pi * 5 * t) * (1 - (-5 * t + 2.5) / 2.5) - 0.3 * np.sin(np.pi * 5 * t) * 2


def Yt2(t):
    return -np.pi * 5 * 0.3 * np.cos(np.pi * 5 * t) * (1 - (5 * t - 2.5) / 2.5) + 0.3 * np.sin(np.pi * 5 * t) * 2


def Yt(t):
    return np.piecewise(t, [t < 0.5, t >= 0.5], [Yt1, Yt2])


class BottomCustomizedRectangle(DomainInputBase):
    """ In our coordinate system, the 'Bottom' actually is the Left side."""
    def __init__(self, h=2, l=None, mapping=(X, Y), Jacobian=(Xt, Yt),):
        """
        The DomainInput describes such a domain:
        ---------------------> y
        |
        |   (0,0)       h
        |     .___________________
        |     |         U        |
        |     \                  |
        |      \                 |
        |      /  L            R |  l
        v     /                  |
        x    /                   |
             \__________D________|
                        h

        The Left side is given by a line: (x, y) = mapping(t) = (X(t), Y(t)), t=[0,1],
        and (X(0), Y(0)) = (0, 0).

        Parameters
        ----------
        h : float
            The height. It must > 0.
        l : float
            The length. It must > 0.
        mapping :
            The mapping.
        Jacobian :
            The Jacobian of the mapping.

        """
        # _____ if l is None, we get is from mapping[0], X _____________________________
        if l is None: l = mapping[0](1)
        # ____ do some checks for the inputs ___________________________________________
        assert isinstance(h, (int, float)) and h > 0, " <BottomCustomizedRectangular> : h must > 0."
        assert isinstance(l, (int, float)) and l > 0, " <BottomCustomizedRectangular> : l must > 0."
        t = np.linspace(0, 1, 20)[1:-1]
        X, Y = mapping
        Xt, Yt = Jacobian
        XND = NumericalDerivative_fx(X, t)
        assert XND.check_derivative(Xt), " <BottomCustomizedRectangular> : Xt is wrong."
        YND = NumericalDerivative_fx(Y, t)
        assert YND.check_derivative(Yt), " <BottomCustomizedRectangular> : Yt is wrong."
        assert np.abs(X(0) - 0) < 1e-8 and np.abs(Y(0) - 0) < 1e-8, " <BottomCustomizedRectangular>."
        assert np.abs(X(1) - l) < 1e-8, " <BottomCustomizedRectangular>."
        # ____ personal properties _____________________________________________________
        self._h_ = h
        self._l_ = l
        self._bottom_mapping_ = mapping
        self._bottom_Jacobian_ = Jacobian
        # _____________ initialize parent ______________________________________________
        super().__init__(domain_name='BottomCustomizedRectangular')
        self.region_corner_coordinates = {'R:R': ((X(0), Y(0)), (X(1), Y(1)), (0, h), (l, h))}
        self.region_edge_types = {'R:R-L': ('customized', self.bottom_mapping, self.bottom_Jacobian), }
        self.boundary_region_edges = {'Upper': "R:R-U", 'Down': "R:R-D", 'Left': "R:R-L", 'Right': "R:R-R"}
        self.region_interpolators = 'transfinite'

        self.region_type_wr2_metric = {'R:R': 'transfinite'}
        self.region_sequence = ('R:R',)

        # ------------------------------------------------------------------------------

    @property
    def h(self):
        return self._h_

    @property
    def l(self):
        return self._l_

    @property
    def bottom_mapping(self):
        return self._bottom_mapping_

    @property
    def bottom_Jacobian(self):
        return self._bottom_Jacobian_



    @classproperty
    def statistic(cls):
        return {'periodic': False,
                'region num': 1,
                'mesh boundary num': 4, # the amount of mesh boundaries (instead of domain boundaries)
                }

    @classproperty
    def random_parameters(cls):
        return {}