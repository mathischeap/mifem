# -*- coding: utf-8 -*-
"""
INTRO

@author: Yi Zhang. Created on Tue May 21 11:57:52 2019
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft
         Delft, Netherlands

"""
import inspect
import numpy as np
from SCREWS.frozen import FrozenOnly
from typing import Dict, Union
from SCREWS.numerical._1d import NumericalDerivative


class DomainInput(FrozenOnly):
    def __init__(self, domain_name='domain without name'):
        self.domain_name = domain_name
        self._ndim_ = 2
        self._region_corner_coordinates_ = None
        self._region_edge_types_ = None
        self._boundary_region_edges_ = None
        self._region_interpolators_ = None

        self._boundary_names_ = None
        self._periodic_boundary_pairs_ = dict()
        self._periodic_boundaries_ = set()
        self._region_sequence_ = None
        self._region_type_wr2_metric_ = None
        self._internal_parameters_ = list()

        INSP = inspect.getfullargspec(self.__init__)
        self.__arg_names___ = INSP[0][1:]
        assert INSP[1] is None and INSP[2] is None, "A domain input class can not have *args and **kwargs."
        assert len(INSP[3]) == len(self.__arg_names___), "A domain input class can only have kwargs."

        self._freeze_self_()

    @property
    def internal_parameters(self):
        """Internal parameters only affect internal metric, does not affect the domain shape."""
        return self._internal_parameters_

    @internal_parameters.setter
    def internal_parameters(self, internal_parameters):
        if isinstance(internal_parameters, list):
            pass
        elif isinstance(internal_parameters, str):
            internal_parameters = [internal_parameters,]
        elif isinstance(internal_parameters, (tuple, set)):
            internal_parameters = list(internal_parameters)
        else:
            raise NotImplementedError(f"internal_parameters = {internal_parameters} not acceptable.")

        assert isinstance(internal_parameters, list), \
            f"please put internal_parameters({internal_parameters}) in list."
        if len(internal_parameters) > 0:
            assert all([ip in self.__arg_names___ for ip in internal_parameters])

        self._internal_parameters_ = internal_parameters

    @property
    def domain_name(self):
        """ Mesh name. """
        return self._domain_name_

    @domain_name.setter
    def domain_name(self, dn):
        assert isinstance(dn, str), " <DomainInput> : domain name needs to be str."
        self._domain_name_ = dn

    @property
    def ndim(self):
        """ dimensions n. """
        return self._ndim_

    @property
    def region_interpolators(self):
        return self._region_interpolators_

    @region_interpolators.setter
    def region_interpolators(self, region_interpolators):
        self._region_interpolators_ = region_interpolators

    def ___PRIVATE_region_name_equirement_checker___(self, regionDict):
        """
        Requirements:
        0). != domain name.
        1). Does not include '-' and '='.
        2). Starts with 'R:'
        """
        for R in regionDict:
            assert R != self.domain_name
            assert '-' not in R and '=' not in R, \
                " <DomainInput> : region name = {} is wrong".format(R)
            assert R[0:2] == 'R:', \
                " <DomainInput> : region name = {} does not start with 'R:'".format(R)


    @property
    def region_corner_coordinates(self):
        """
        Store the regions 4 corners' coordinates.

        Returns
        -------
        region_coordinates : dict
            A dict whose keys represent the region names, and values represent
            the coordinates of region corner points.

            In 2D: (UL, DL, UR, DR).

            L: Left, R: Right, U: Upper, D: Down

        """
        return self._region_corner_coordinates_

    @region_corner_coordinates.setter
    def region_corner_coordinates(self, _dict_):
        assert isinstance(_dict_, dict), " <DomainInput> : region_coordinates needs to be a dict."

        self.___PRIVATE_region_name_equirement_checker___(_dict_)

        for R in _dict_:
            assert isinstance(R, str) and R != '' and '-' not in R, " <DomainInput> : region name = {} is wrong".format(
                R)
            assert np.shape(_dict_[R])[0] == 4, \
                " <DomainInput> : region_coordinates[{}]={} is wrong.".format(R, _dict_[R])
        self._region_corner_coordinates_ = _dict_

    @property
    def region_edge_types(self):
        """
        Store the regions's boundaries' types.

        Returns
        -------
        region_boundary_type : dict
            A dict that contains the regions boundary info. The keys indicate
            the region boundary, the value indicate the info. value[0] indicate
            the type, value[1:] indicate the rest info which will be parsed
            into full information. The not mentioned region boundaries will be
            set into default type: ('plane',)

            Notice that the value will be send to edge_geometry eventaully. And
            if these is info (value[1:]) to be parsed, it will be done there in
            edge_geometry. And the value is stored in the
            `edge_geometry.edge_types`.

        """
        return self._region_edge_types_

    @region_edge_types.setter
    def region_edge_types(self, _dict_):
        assert self.region_corner_coordinates is not None, " <DomainInput> : please first set region_coordinates."
        assert isinstance(_dict_, dict), " <DomainInput> : region_boundary_type needs to be a dict."
        for item in _dict_:
            R, S = item.split('-')
            assert R in self.region_corner_coordinates and S in ('U', 'D', 'L', 'R'), \
                " <DomainInput> : region_edge_type key {} is wrong.".format(item)
        self._region_edge_types_ = _dict_



    def ___PRIVATE_boundary_name_requirement_checker___(self, boundaryRegionSidesDict):
        """
        Requirements:
        0). != domain name.
        1). Is String and is not empty. Does not include '-' and '='.
        2). Can not start with 'R:' (So it must be different from region names).
        3). Length > 1
        4). Does not contain any number
        """
        for boundary_name in boundaryRegionSidesDict.keys():
            assert boundary_name != self.domain_name
            assert isinstance(boundary_name, str) and boundary_name != '' \
                   and '-' not in boundary_name and '=' not in boundary_name, \
                " <DomainInput> : boundary_name = {} is wrong.".format(boundary_name)
            assert boundary_name[0:2] not in ('R:',), \
                " <DomainInput> : boundary_name = {} wrong.".format(boundary_name)
            assert len(boundary_name) > 1, \
                " <DomainInput> : boundary_name = {} is too short.".format(boundary_name)
            for i in '0123456789':
                assert i not in boundary_name, \
                    " <DomainInput> : boundary_name = {} has number {}.".format(boundary_name, i)


    @property
    def boundary_region_edges(self):
        """
        Store the domain boundary infomation.

        Returns
        -------
        domain_boundary : dict
            For example:
                {'Down': ("Body_center-D", 'Body_back-D', ...),
                 'West': ("Body_center-R", 'Body_back-R', ...),
                 ......}
            This means we have domain boundaries: South, West and so on.

        """
        return self._boundary_region_edges_

    @boundary_region_edges.setter
    def boundary_region_edges(self, _dict_):
        assert self.region_corner_coordinates is not None, " <DomainInput> : please first set region_coordinates."
        assert isinstance(_dict_, dict), " <DomainInput> : domain_boundary needs to be a dict."

        self.___PRIVATE_boundary_name_requirement_checker___(_dict_)

        for boundary_names in _dict_.keys():
            assert isinstance(boundary_names, str) and boundary_names != '' and '-' not in boundary_names, \
                " <DomainInput> : boundary_names = {} is wrong.".format(boundary_names)
            assert boundary_names not in self.region_corner_coordinates.keys(), \
                " <DomainInput>: boundary_names={} is taken by one of the regions.".format(boundary_names)
        for item in _dict_:
            if isinstance(_dict_[item], str):
                _dict_[item] = (_dict_[item],)
            if isinstance(_dict_[item], list) or isinstance(_dict_[item], tuple):
                for item_i in _dict_[item]:
                    R, S = item_i.split('-')
                    assert R in self.region_corner_coordinates and S in ('U', 'D', 'L', 'R'), \
                        " <DomainInput> : domain_boundary[{}]={} is wrong.".format(item, _dict_[item])
            else:
                raise Exception(" <DomainInput> : boundary_region_edges input value accepts only str, tuple of list.")
        self._boundary_region_edges_ = _dict_
        self._boundary_names_ = list(_dict_.keys())






    def ___PRIVATE_periodic_boundary_requirement_checker___(self, pBd):
        """
        Here we only do a simple check. We make sure that the keys are in format of:
        0). boundary_name_1=boundary_name_2.
        1). A boundary name at most appear in one pair.

        """
        assert isinstance(pBd, dict)
        bnPOOL = set()
        for pair in pBd:
            assert '=' in pair
            bn1, bn2 = pair.split('=')
            lengthPOOL = len(bnPOOL)
            assert bn1 in self._boundary_names_ and bn2 in self._boundary_names_
            bnPOOL.add(bn1)
            bnPOOL.add(bn2)
            newLengthPOOL = len(bnPOOL)
            assert newLengthPOOL == lengthPOOL + 2, "Boundary(s) used for multiple periodic pairs!"
        self._periodic_boundaries_ = bnPOOL


    @property
    def periodic_boundary_pairs(self):
        return self._periodic_boundary_pairs_

    @periodic_boundary_pairs.setter
    def periodic_boundary_pairs(self, pBd):
        """ """
        self.___PRIVATE_periodic_boundary_requirement_checker___(pBd)
        self._periodic_boundary_pairs_ = pBd

    @property
    def periodic_boundaries(self):
        """(set) Return a set of all boundary names those involved in the periodic boundary setting."""
        return self._periodic_boundaries_

    @property
    def periodic_boundaries_involved_regions(self):
        """"""
        regions = set()
        for pb in self.periodic_boundaries:
            region_sides = self.boundary_region_edges[pb]
            for rs in region_sides:
                rn = rs.split('-')[0]
                if rn not in regions:
                    regions.add(rn)
        return regions





    @property
    def region_sequence(self):
        """
        This will fix the sequence of regions by fix their names in property
        region_names or regions.names. This is very important for numbering. Sometimes, a bad
        region sequence can make the numbering wrong.

        """
        return self._region_sequence_

    @region_sequence.setter
    def region_sequence(self, rS: tuple):
        assert len(rS) == len(self.region_corner_coordinates.keys())
        assert all([rSi in self.region_corner_coordinates for rSi in rS]) & \
            all([rSi in rS for rSi in self.region_corner_coordinates.keys()]), \
            f"region_sequence={rS} has invalid region name(s)."
        self._region_sequence_ = rS

    @property
    def region_type_wr2_metric(self):
        return self._region_type_wr2_metric_

    @region_type_wr2_metric.setter
    def region_type_wr2_metric(self, rTwr2M: Union[str, Dict[str, str]]):
        if isinstance(rTwr2M, str):
            _D_ = dict()
            for region_name in self.region_corner_coordinates:
                _D_[region_name] = rTwr2M
            rTwr2M = _D_
        assert isinstance(rTwr2M, dict), "region_type_wr2_metric needs to be a dictionary."
        for key in rTwr2M:
            assert key in self.region_corner_coordinates, f"Region name={key} not valid."
        self._region_type_wr2_metric_ = rTwr2M



# ------------- particular domain inputs -------------------------------------------------------------




class Crazy(DomainInput):
    def __init__(self, c=0, bounds=((0, 1), (0, 1))):
        assert np.shape(bounds)[0] == 2, " <Crazy> : bounds dimension is wrong."
        for i in range(2):
            assert np.shape(bounds[i]) == (2,) and bounds[i][1] > bounds[i][0], \
                " <CrazyDomain> : bounds[{}]=={} is wrong.".format(i, bounds[i])
        self._bounds_ = bounds
        self._c_ = c
        super().__init__(domain_name='Crazy')
        x0, x1 = bounds[0]
        y0, y1 = bounds[1]
        # _____________ standard inputs ________________________________________________
        self.region_corner_coordinates = {'R:R': ((x0, y0), (x1, y0), (x0, y1), (x1, y1))}
        self.region_edge_types = {'R:R-L': ('free',)}
        self.boundary_region_edges = {'Upper': "R:R-U", 'Down': "R:R-D", 'Left': "R:R-L", 'Right': "R:R-R"}
        self.region_interpolators = 'crazy'

        self.region_type_wr2_metric = {'R:R': 'crazy'}
        self.region_sequence = ('R:R',)
        self.internal_parameters = ['c', ]  # has to be defined after the super().__init__
        # ------------------------------------------------------------------------------

    @property
    def bounds(self):
        return self._bounds_

    @property
    def c(self):
        return self._c_



class CrazyPeriodic(DomainInput):
    def __init__(self, c=0, bounds=((0, 1), (0, 1))):
        assert np.shape(bounds)[0] == 2, " <Crazy_Periodic> : bounds dimension is wrong."
        for i in range(2):
            assert np.shape(bounds[i]) == (2,) and bounds[i][1] > bounds[i][0], \
                " <Crazy_Periodic> : bounds[{}]=={} is wrong.".format(i, bounds[i])
        self._bounds_ = bounds
        self._c_ = c
        super().__init__(domain_name='CrazyPeriodic')
        x0, x1 = bounds[0]
        y0, y1 = bounds[1]
        # _____________ standard inputs ________________________________________________
        self.region_corner_coordinates = {'R:R': ((x0, y0), (x1, y0), (x0, y1), (x1, y1))}
        self.region_edge_types = {'R:R-L': ('free',)}
        self.boundary_region_edges = {'Upper': "R:R-U", 'Down': "R:R-D", 'Left': "R:R-L", 'Right': "R:R-R"}
        self.region_interpolators = 'crazy'

        self.periodic_boundary_pairs = {'Upper=Down': 'regular',
                                        'Left=Right'  : 'regular',}
        self.region_type_wr2_metric = {'R:R': 'crazy'}
        self.region_sequence = ('R:R',)
        self.internal_parameters = ['c', ]  # has to be defined after the super().__init__
        # ------------------------------------------------------------------------------

    @property
    def bounds(self):
        return self._bounds_

    @property
    def c(self):
        return self._c_




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


class BottomCustomizedRectangle(DomainInput):
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
        XND = NumericalDerivative(X, t)
        assert XND.check_derivative(Xt), " <BottomCustomizedRectangular> : Xt is wrong."
        YND = NumericalDerivative(Y, t)
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



class Quadrangle(DomainInput):
    """ A quadrangle of four straight edges"""
    def __init__(self, p_UL=(0,0), p_DL=(1,0), p_UR=(1,1), p_DR=(2,1)):
        super().__init__(domain_name='Quadrangle')
        self.region_corner_coordinates = {'R:R': (p_UL, p_DL, p_UR, p_DR)}
        self.region_edge_types = {}
        self.boundary_region_edges = {'Upper': "R:R-U", 'Down': "R:R-D", 'Left': "R:R-L", 'Right': "R:R-R"}
        self.region_interpolators = 'transfinite'
        self.region_type_wr2_metric = {'R:R': 'transfinite'}
        self.region_sequence = ('R:R',)





class CylinderInChannel(DomainInput):
    """ Just like the class name say, this is a cylinder in channel domain."""
    def __init__(self, h=1.5, r=0.5, li=0.75, lo=2.25):
        """
        The DomainInput describes such a domain:
                            y
                            ^
                            |
         ____li_____._______|________.___________lo___________
        |           .       |        .                        |
        |           .     __|__      .                        |
        |           .    /  |  \     .                        |
        |h          .h  |   .-r-|----.------------------------|-----> x
        |           .    \_____/     .                        |
        |           .                .                        |
        |____li_____._______h________.___________lo___________|

        The domian is divided into following regions:
         ___________._______________.________________________
        |           .\      Ru     /.                        |
        |           . \   _____   / .                        |
        |           .  \ /     \ /  .                        |
        |     Ri    .Rl |   .-r-| Rr.           Ro           |
        |           .  / \_____/ \  .                        |
        |           . /     Rd    \ .                        |
        |___________./_____________\.________________________|


        Parameters
        ----------
        h : float
            The height. It must > 0.
        r : float
            The radius of the cylinder. It must > 0.
        li :
            The Inlet Length. It must > 0.
        lo :
            The Outlet Length. It must > 0.

        """
        # ____ do some checks first ____________________________________________________
        assert isinstance(h, (int, float)) and h > 0, " <CylinderInChannel> : need h > 0."
        assert isinstance(r, (int, float)) and r > 0, " <CylinderInChannel> : need r > 0."
        assert isinstance(li, (int, float)) and r > 0, " <CylinderInChannel> : need li > 0."
        assert isinstance(lo, (int, float)) and r > 0, " <CylinderInChannel> : need lo > 0."
        assert r < h / 2, " <CylinderInChannel> : need r < h/2."
        # ____ personal properties _____________________________________________________
        self._h_ = h
        self._r_ = r
        self._li_ = li
        self._lo_ = lo
        # ______ initialiing parent ____________________________________________________
        hsqr_r = 0.5 * np.sqrt(r)
        super().__init__(domain_name='CylinderInChannel')
        self.region_corner_coordinates = {
            'R:Ru': ((-hsqr_r, hsqr_r), (hsqr_r, hsqr_r), (-h / 2, h / 2), (h / 2, h / 2)),
            'R:Rr': ((hsqr_r, hsqr_r), (hsqr_r, -hsqr_r), (h / 2, h / 2), (h / 2, -h / 2)),
            'R:Rd': ((hsqr_r, -hsqr_r), (-hsqr_r, -hsqr_r), (h / 2, -h / 2), (-h / 2, -h / 2)),
            'R:Rl': ((-hsqr_r, -hsqr_r), (-hsqr_r, hsqr_r), (-h / 2, -h / 2), (-h / 2, h / 2)),
            'R:Ri': ((-h / 2, -h / 2), (-h / 2, h / 2), (-h / 2 - li, -h / 2), (-h / 2 - li, h / 2)),
            'R:Ro': ((h / 2, h / 2), (h / 2, -h / 2), (h / 2 + lo, h / 2), (h / 2 + lo, -h / 2))}
        self.region_edge_types = {'R:Ru-L': ('acw', (0, 0)),
                                  'R:Rl-L': ('acw', (0, 0)),
                                  'R:Rd-L': ('acw', (0, 0)),
                                  'R:Rr-L': ('acw', (0, 0)), }
        self.boundary_region_edges = {'Upper': ('R:Ri-D', 'R:Ru-R', 'R:Ro-U'),
                                      'Down': ('R:Ri-U', 'R:Rd-R', 'R:Ro-D'),
                                      'Left': 'R:Ri-R',
                                      'Right': 'R:Ro-R',
                                      'Internal': ('R:Ru-L', 'R:Rl-L', 'R:Rd-L', 'R:Rr-L')}
        self.region_interpolators = 'transfinite'

        self.region_type_wr2_metric = 'transfinite'
        # self.region_sequence = ('R:R',)
        # self._internal_parameters_ = None

        # ------------------------------------------------------------------------------

    @property
    def h(self):
        return self._h_

    @property
    def r(self):
        return self._r_

    @property
    def li(self):
        return self._li_

    @property
    def lo(self):
        return self._lo_


class CircleHolePlate1(DomainInput):
    """ """
    def __init__(self, hx=1.5, hy=None, r=0.5):
        """
        The DomainInput describes such a domain:
                y
                ^
           hx/2 |  hx/2
         _______|________
        |       |        |
        |     __|__      |hy/2
        |    /  |  \     |
        |hy |   .-r-|----|---> x
        |    \_____/     |
        |                |hy/2
        |_______hx_______|

        The domian is divided into following regions:
         _______________
        |\      Ru     /|
        | \   _____   / |
        |  \ /     \ /  |
        |Rl |   .-r-| Rr|
        |  / \_____/ \  |
        | /     Rd    \ |
        |/_____________\|


        Parameters
        ----------

        """
        # ____ parse inputs ____________________________________________________________
        if hy is None: hy = hx
        # ____ checks __________________________________________________________________
        assert hx > 0 and hy > 0, " <HolePlate> : hx={}, hy={} illegal.".format(hx, hy)
        assert r < np.min((hx, hy)), " <HolePlate> : r={} too large.".format(r)
        # ______ initialiing parent ____________________________________________________
        hsqr_r = 0.5 * np.sqrt(r)
        super().__init__(domain_name='CircleHolePlate1')
        self.region_corner_coordinates = {
            'R:Ru': ((-hsqr_r, hsqr_r), (hsqr_r, hsqr_r), (-hx / 2, hy / 2), (hx / 2, hy / 2)),
            'R:Rr': ((hsqr_r, hsqr_r), (hsqr_r, -hsqr_r), (hx / 2, hy / 2), (hx / 2, -hy / 2)),
            'R:Rd': ((hsqr_r, -hsqr_r), (-hsqr_r, -hsqr_r), (hx / 2, -hy / 2), (-hx / 2, -hy / 2)),
            'R:Rl': ((-hsqr_r, -hsqr_r), (-hsqr_r, hsqr_r), (-hx / 2, -hy / 2), (-hx / 2, hy / 2))}
        self.region_edge_types = {'R:Ru-L': ('acw', (0, 0)),
                                  'R:Rl-L': ('acw', (0, 0)),
                                  'R:Rd-L': ('acw', (0, 0)),
                                  'R:Rr-L': ('acw', (0, 0)), }
        self.boundary_region_edges = {'Upper': 'R:Ru-R',
                                      'Down': 'R:Rd-R',
                                      'Left': 'R:Rl-R',
                                      'Right': 'R:Rr-R',
                                      'Internal': ('R:Ru-L', 'R:Rl-L', 'R:Rd-L', 'R:Rr-L')}
        self.region_interpolators = 'transfinite'

        self.region_type_wr2_metric = 'transfinite'
        self.internal_parameters = list()
        # ------------------------------------------------------------------------------


class CircleHolePlate2(DomainInput):
    """ """
    def __init__(self, hx=2, hy=None, r=0.5):
        """
        The DomainInput describes such a domain:

              ________hy_________
             |                   |
             |        ___        |hx/2
             |       /   \       |
           hx|      |  .r |      |-----> y
             |       \___/       |
             |                   |
             |___________________|
                 hy/2  |
                       |
                       |
                       v x

        The domain is divided into following regions:
                       U
              ___________________
             |       |R_U|       |
             | R_UL  |___|  R_UR |
             |-------/   \-------|
           L | R_L  |  .  | R_R  | R
             |-------\___/-------|
             | R_DL  |R_D|  R_DR |
             |_______|___|_______|
                       D

        The center of the circle hole is at (0, 0).

        Parameters
        ----------

        """
        # ____ parse inputs ____________________________________________________________
        if hy is None: hy = hx
        # ____ checks __________________________________________________________________
        assert hx > 0 and hy > 0, " <HolePlate> : hx={}, hy={} illegal.".format(hx, hy)
        assert r < np.min((hx, hy)), " <HolePlate> : r={} too large.".format(r)
        # _____________ standard inputs ________________________________________________
        super().__init__(domain_name='CircleHolePlate2')
        sr = np.sqrt(2) * r / 2
        self.region_corner_coordinates = {
            'R:R_UL': ((-hx / 2, -hy / 2), (-sr, -hy / 2), (-hx / 2, -sr), (-sr, -sr)),
            'R:R_L': ((-sr, -hy / 2), (sr, -hy / 2), (-sr, -sr), (sr, -sr)),
            'R:R_DL': ((sr, -hy / 2), (hx / 2, -hy / 2), (sr, -sr), (hx / 2, -sr)),
            'R:R_U': ((-hx / 2, -sr), (-sr, -sr), (-hx / 2, sr), (-sr, sr)),
            'R:R_D': ((sr, -sr), (hx / 2, -sr), (sr, sr), (hx / 2, sr)),
            'R:R_UR': ((-hx / 2, sr), (-sr, sr), (-hx / 2, hy / 2), (-sr, hy / 2)),
            'R:R_R': ((-sr, sr), (sr, sr), (-sr, hy / 2), (sr, hy / 2)),
            'R:R_DR': ((sr, sr), (hx / 2, sr), (sr, hy / 2), (hx / 2, hy / 2))}
        self.region_edge_types = {'R:R_L-R': ('aacw', (0, 0)),
                                  'R:R_U-D': ('acw', (0, 0)),
                                  'R:R_D-U': ('aacw', (0, 0)),
                                  'R:R_R-L': ('acw', (0, 0)), }
        self.boundary_region_edges = {'Upper': ("R:R_UL-U", 'R:R_U-U', "R:R_UR-U"),
                                      'Down': ("R:R_DL-D", 'R:R_D-D', "R:R_DR-D"),
                                      'Left': ("R:R_UL-L", 'R:R_L-L', "R:R_DL-L"),
                                      'Right': ("R:R_UR-R", 'R:R_R-R', "R:R_DR-R"),
                                      'Internal': ("R:R_L-R", "R:R_U-D", "R:R_D-U", 'R:R_R-L')}
        self.region_interpolators = 'transfinite'

        self.region_type_wr2_metric = 'transfinite'
        self.internal_parameters = list()

if __name__ == "__main__":
    C = Crazy(c=0)