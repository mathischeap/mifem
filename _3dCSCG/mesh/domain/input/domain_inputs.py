# -*- coding: utf-8 -*-
"""
Hereby, we define a class for defining a domain. Each instance of such a class
define stores the input information for a domain, 2D or 3D. The mesh_generator
will take such a instance as input and then call the correct generator to 
generate a based on the domain.

@author: Yi Zhang, Created on Mon Sep  3 12:13:28 2018
         Aerodynamics, AE
         TU Delft
"""
import inspect
import numpy as np
from typing import Dict
from SCREWS.frozen import FrozenOnly
from SCREWS.miscellaneous import check_no_splcharacter

class _3dDomainInput(FrozenOnly):
    def __init__(self, domain_name='domain without name'):
        self.domain_name = domain_name
        self._ndim_ = 3
        self._region_corner_coordinates_ = None
        self._region_side_types_ = None
        self._boundary_region_sides_ = None
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
    def region_interpolators(self, rip):
        assert isinstance(rip, (dict, str)), "region_interpolators needs to be str or dict."
        self._region_interpolators_ = rip


    def ___PRIVATE_region_name_requirement_checker___(self, regionDict):
        """
        Requirements:
        0). != domain name.
        1). Does not include '-' and '='.
        2). Starts with 'R:'
        """
        for R in regionDict:
            assert R != self.domain_name
            assert check_no_splcharacter(R[2:]), f"region name {R} has special characters."
            assert '-' not in R and '=' not in R, \
                " <DomainInput> : region name = {} is wrong".format(R)
            assert R[0:2] == 'R:', \
                " <DomainInput> : region name = {} does not start with 'R:'".format(R)


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
            assert check_no_splcharacter(boundary_name), f"boundary name {boundary_name} has special characters."
            assert isinstance(boundary_name, str) and boundary_name != '' \
                   and '-' not in boundary_name and '=' not in boundary_name \
                   and '|' not in boundary_name, \
                " <DomainInput> : boundary_name = {} is wrong.".format(boundary_name)
            assert boundary_name[0:2] not in ('R:',), \
                " <DomainInput> : boundary_name = {} wrong.".format(boundary_name)
            assert len(boundary_name) > 1, \
                " <DomainInput> : boundary_name = {} is too short.".format(boundary_name)
            for i in '0123456789':
                assert i not in boundary_name, \
                    " <DomainInput> : boundary_name = {} has number {}.".format(boundary_name, i)
        self._boundary_names_ = list(boundaryRegionSidesDict.keys())





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
            region_sides = self.boundary_region_sides[pb]
            for rs in region_sides:
                rn = rs.split('-')[0]
                if rn not in regions:
                    regions.add(rn)
        return regions


    @property
    def region_corner_coordinates(self):
        """
        Store the regions 8 corners' coordinates.
        
        Returns
        -------
        region_coordinates : dict
            A dict whose keys represent the region names, and values represent
            the coordinates of region corner points.
            
            In 3D: (NWB, SWB, NEB, SEB, NWF, SWF, NEF, SEF).
            S: South, N: North, W: West, E:East, B: Back, F: Front.
            
        """
        return self._region_corner_coordinates_
    
    @region_corner_coordinates.setter
    def region_corner_coordinates(self, _dict_):
        assert isinstance(_dict_, dict), " <DomainInput> : region_coordinates needs to be a dict."
        self.___PRIVATE_region_name_requirement_checker___(_dict_)
        for R in _dict_:
            assert np.shape(_dict_[R])[0] == 8, \
                " <DomainInput> : region_coordinates[{}]={} is wrong.".format(R, _dict_[R])
        self._region_corner_coordinates_ = _dict_

    @property
    def region_side_types(self):
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
            
            Notice that the value will be send to side_geometry eventually. And
            if these is info (value[1:]) to be parsed, it will be done there in
            side_geometry. And the value is stored in the 
            `side_geometry.side_types`.
            
        """
        return self._region_side_types_
    
    @region_side_types.setter
    def region_side_types(self, _dict_):
        assert self.region_corner_coordinates is not None, \
            " <DomainInput> : please first set region_coordinates."
        assert isinstance(_dict_, dict), \
            " <DomainInput> : region_boundary_type needs to be a dict."
        for item in _dict_:
            R, S = item.split('-')
            assert R in self.region_corner_coordinates and S in ('N', 'S', 'W', 'E', 'B', 'F'), \
                " <DomainInput> : region_side_type key {} is wrong.".format(item)
        self._region_side_types_ = _dict_

    @property
    def boundary_region_sides(self):
        """
        Store the domain boundary information.
        
        Returns
        -------
        domain_boundary : dict
            For example:
                {'South': ("Body_center-S", 'Body_back-S', "Body_front-S"), 
                 'West': ("Body_center-W", 'Body_back-W', "Body_front-W", "Body_north-W"),
                 ......}
            This means we have three domain boundaries: South, West and so on.
            The South contains the South side of region: Body_center, the South
            side of region: Body_back and the South side of region: Body_front.
            
        """
        return self._boundary_region_sides_
    
    @boundary_region_sides.setter
    def boundary_region_sides(self, _dict_):
        assert self.region_corner_coordinates is not None, \
            " <DomainInput> : please first set region_coordinates."
        assert isinstance(_dict_, dict), " <DomainInput> : domain_boundary needs to be a dict."

        self.___PRIVATE_boundary_name_requirement_checker___(_dict_)

        for item in _dict_:
            if isinstance(_dict_[item], str):
                _dict_[item] = (_dict_[item],)
            regionSidePool = set()
            if isinstance(_dict_[item], list) or isinstance(_dict_[item], tuple):
                assert len(_dict_[item]) > 0, "<DomainInput> : boundary {} is empty.".format(item)
                for item_i in _dict_[item]:
                    assert item_i not in regionSidePool, f"regionSide:{item_i} appears in multiple boundaries!"
                    regionSidePool.add(item_i)
                    R, S = item_i.split('-')
                    assert R in self.region_corner_coordinates and S in ('N', 'S', 'W', 'E', 'B', 'F'), \
                        " <DomainInput> : domain_boundary[{}]={} is wrong.".format(item, _dict_[item])
            else:
                raise Exception(" <DomainInput> : boundary_region_sides input value accepts" + \
                                " only str, tuple of list.")
        self._boundary_region_sides_ = _dict_

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
    def region_type_wr2_metric(self, rTwr2M: Dict[str, str]):
        assert isinstance(rTwr2M, dict), "region_type_wr2_metric needs to be a dictionary."
        for key in rTwr2M:
            assert key in self.region_corner_coordinates, f"Region name={key} not valid."
        self._region_type_wr2_metric_ = rTwr2M




# ----- particular domains below ---------------------------------------------------------------------------------------



class BridgeArchCracked(_3dDomainInput):
    def __init__(self, domain_name="BridgeArch",
        l=4, w=1, h=1.5, r=5 / 2, d=0.25,
        region_interpolators='bridge_arch_cracked'):
        """
        Parameters
        ----------
        r : float
            Arch radius.
        d : float
            Uncracked depth.

        """
        self._l_, self._w_, self._h_, self._r_ = l, w, h, r
        assert r ** 2 - l ** 2 / 4 > 0, " <BridgeArch> : parameters do not fit an arch."
        alpha = np.sqrt(r ** 2 - l ** 2 / 4)
        beta = h + alpha - r
        assert beta > 0, " <BridgeArch> : bridge has zero thickness."
        self._center_ = (h + alpha, l / 2)
        self._alpha_ = alpha
        self._beta_ = beta
        self._d_ = d
        assert d < beta, " <BridgeArchCracked> : Bridge already breaks into two."
        A = (0, 0, 0)
        B = (h, 0, 0)
        C = (0, l / 2, 0)
        Dl = (beta, l / 2, 0)
        Dr = (beta, l / 2, 0)
        E = (0, 0, w)
        F = (h, 0, w)
        G = (0, l / 2, w)
        Hl = (beta, l / 2, w)
        Hr = (beta, l / 2, w)
        I = (0, l, 0)
        J = (h, l, 0)
        K = (0, l, w)
        L = (h, l, w)
        rd = beta - d
        M = (rd * beta / h, 0, 0)
        N = (rd * beta / h, 0, w)
        O = (rd, l / 2, 0)
        P = (rd, l / 2, w)
        Q = (rd * beta / h, l, 0)
        R = (rd * beta / h, l, w)
        super().__init__(domain_name=domain_name)
        self.region_corner_coordinates = {'R:R_left_up': (A, M, C, O, E, N, G, P),
                                          'R:R_left_down': (M, B, O, Dl, N, F, P, Hl),
                                          'R:R_right_up': (C, O, I, Q, G, P, K, R),
                                          'R:R_right_down': (O, Dr, Q, J, P, Hr, R, L)}
        self.region_side_types = {} # all plane
        self.boundary_region_sides = {
            'Left_Floor': ("R:R_left_up-N",),
            'Right_Floor': ('R:R_right_up-N',),
            'Back': ("R:R_left_up-B", "R:R_left_down-B", 'R:R_right_up-B', 'R:R_right_down-B',),
            'Front': ("R:R_left_up-F", "R:R_left_down-F", 'R:R_right_up-F', 'R:R_right_down-F',),
            'Bottom': ("R:R_left_down-S", 'R:R_right_down-S',),
            'Left_Wall': ('R:R_left_up-W', 'R:R_left_down-W',),
            'Right_Wall': ('R:R_right_up-E', 'R:R_right_down-E',),
            'Crack': ('R:R_left_down-E', 'R:R_right_down-W')
        }
        self.region_interpolators = region_interpolators

    @property
    def l(self):
        return self._l_

    @property
    def w(self):
        return self._w_

    @property
    def h(self):
        return self._h_

    @property
    def r(self):
        """Arch radius."""
        return self._r_

    @property
    def d(self):
        """Crack depth."""
        return self._d_

    @property
    def center(self):
        return self._center_

    @property
    def beta(self):
        return self._beta_


class Crazy(_3dDomainInput):
    """A "crazy" 3d rectangular domain's input class whose inside is distorted with the "crazy" mapping."""

    def __init__(self, c=0, bounds=((0, 1), (0, 1), (0, 1)), domain_name="Crazy3D"):
        assert np.shape(bounds)[0] == 3, " <CrazyDomain> : bounds dimension is wrong."
        for i in range(3):
            assert np.shape(bounds[i]) == (2,) and bounds[i][1] > bounds[i][0], \
                " <CrazyDomain> : bounds[{}]=={} is wrong.".format(i, bounds[i])
        self._bounds_ = bounds
        self._c_ = c
        super().__init__(domain_name=domain_name)

        x0, x1 = bounds[0]
        y0, y1 = bounds[1]
        z0, z1 = bounds[2]
        assert x1 > x0
        assert y1 > y0
        assert z1 > z0
        self.region_corner_coordinates = {'R:R': ((x0, y0, z0), (x1, y0, z0), (x0, y1, z0), (x1, y1, z0),
                                                  (x0, y0, z1), (x1, y0, z1), (x0, y1, z1), (x1, y1, z1))}
        self.region_side_types = {'R:R-S': ('free',)}
        self.boundary_region_sides = {'South': "R:R-S", 'North': "R:R-N",
                                      'West': "R:R-W", 'East': "R:R-E",
                                      'Back': "R:R-B", 'Front': "R:R-F"}
        self.region_interpolators = 'crazy'
        self.region_type_wr2_metric = {'R:R': 'crazy'}
        self.internal_parameters = ['c', ]  # has to be defined after the super().__init__

    @property
    def bounds(self):
        return self._bounds_

    @property
    def c(self):
        return self._c_


class CrazyPeriodic(_3dDomainInput):
    """A "crazy" 3d rectangular domain's input class whose inside is distorted with the "crazy" mapping.
    """
    def __init__(self, c=0, bounds=((0, 1), (0, 1), (0, 1)), domain_name="CrazyPeriodic3D"):
        assert np.shape(bounds)[0] == 3, " <CrazyDomain> : bounds dimension is wrong."
        for i in range(3):
            assert np.shape(bounds[i]) == (2,) and bounds[i][1] > bounds[i][0], \
                " <CrazyDomain> : bounds[{}]=={} is wrong.".format(i, bounds[i])
        self._bounds_ = bounds
        self._c_ = c
        super().__init__(domain_name=domain_name)

        x0, x1 = bounds[0]
        y0, y1 = bounds[1]
        z0, z1 = bounds[2]
        assert x1 > x0 and y1 > y0 and z1 > z0
        self.region_corner_coordinates = {'R:R': ((x0, y0, z0), (x1, y0, z0), (x0, y1, z0), (x1, y1, z0),
                                                  (x0, y0, z1), (x1, y0, z1), (x0, y1, z1), (x1, y1, z1))}
        self.region_side_types = {'R:R-S': ('free',)}
        self.boundary_region_sides = {'North': "R:R-N", 'South': "R:R-S",
                                      'West': "R:R-W" , 'East': "R:R-E" ,
                                      'Back': "R:R-B" , 'Front': "R:R-F"}
        self.region_interpolators = {'R:R':'crazy'}
        self.periodic_boundary_pairs = {'South=North': 'regular',
                                        'West=East'  : 'regular',
                                        'Back=Front' : 'regular'}
        self.region_sequence = ('R:R',)
        self.region_type_wr2_metric = {'R:R': 'crazy'}
        self.internal_parameters = ['c',] # has to be defined after the super().__init__

    @property
    def bounds(self):
        return self._bounds_

    @property
    def c(self):
        return self._c_



class Periodic_Square_Channel(_3dDomainInput):
    """A periodic square channel domain.

    The domain is periodic along x-direction. The other four sides are normal boundaries.

    ^ y
    |                           l
    |_____________________________________________________________
    |                                                            |
    |                                                            |
    |w                                                           |w
    |------------------------------------------------------------|----------> x
    |                                                            |
    |                                                            |
    |                                                            |
    ——————————————————————————————————————————————————————————————
                                l
    """

    def __init__(self, l=2, w=1, h=1, domain_name="Periodic-Square-Channel"):
        assert l > 0 and w > 0 and h > 0, f"l, w, h = {l}, {w}, {h} is wrong."

        self._lwh_ = [l, w, h]

        super().__init__(domain_name=domain_name)

        x0 = 0
        x1 = l
        y0 = - w / 2
        y1 = + w / 2
        z0 = - h / 2
        z1 = + h / 2

        self.region_corner_coordinates = {'R:R': ((x0, y0, z0), (x1, y0, z0), (x0, y1, z0), (x1, y1, z0),
                                                  (x0, y0, z1), (x1, y0, z1), (x0, y1, z1), (x1, y1, z1))}
        self.region_side_types = dict() # all plane
        self.boundary_region_sides = {'wXm': "R:R-N", 'wXp': "R:R-S",
                                      'wYm': "R:R-W", 'wYp': "R:R-E",
                                      'wZm': "R:R-B", 'wZp': "R:R-F"}
        self.region_interpolators = {'R:R': 'transfinite'}
        self.periodic_boundary_pairs = {'wXm=wXp': 'regular',}
        self.region_type_wr2_metric = {'R:R': 'transfinite'}
        self.internal_parameters = list()  # has to be defined after the super().__init__

    @property
    def lwh(self):
        return self._lwh_



class Parallel_Wall_Channel(_3dDomainInput):
    """A channel of two parallel walls.

    Stream-wise: x-direction (length).
    Span-wise: y-direction (width).
    height: z.

    The domain is periodic along x-direction. The other four sides are normal boundaries.

    ^ z
    |                           l
    |_____________________________________________________________
    |                                                            |
    |                                                            |
    |h                                                           |h
    |------------------------------------------------------------|----------> x
    |                                                            |
    |                                                            |
    |                                                            |
    ——————————————————————————————————————————————————————————————
                                l
    """

    def __init__(self, l=2, w=1, h=1, domain_name="Parallel-Wall-Channel"):
        assert l > 0 and w > 0 and h > 0, f"l, w, h = {l}, {w}, {h} is wrong."

        self._lwh_ = [l, w, h]

        super().__init__(domain_name=domain_name)

        x0 = 0
        x1 = l
        y0 = - w / 2
        y1 = + w / 2
        z0 = - h / 2
        z1 = + h / 2

        self.region_corner_coordinates = {'R:R': ((x0, y0, z0), (x1, y0, z0), (x0, y1, z0), (x1, y1, z0),
                                                  (x0, y0, z1), (x1, y0, z1), (x0, y1, z1), (x1, y1, z1))}
        self.region_side_types = dict() # all plane
        self.boundary_region_sides = {'Inlet': "R:R-N", 'Outlet': "R:R-S",
                                      'SpanM': "R:R-W", 'SpanP': "R:R-E",
                                      'wZm': "R:R-B", 'wZp': "R:R-F"}
        self.region_interpolators = {'R:R': 'transfinite'}
        self.periodic_boundary_pairs = {'Inlet=Outlet': 'regular',
                                        'SpanM=SpanP': 'regular',}
        self.region_type_wr2_metric = {'R:R': 'transfinite'}
        self.internal_parameters = list()  # has to be defined after the super().__init__

    @property
    def lwh(self):
        return self._lwh_



class Lid_Driven_Cavity(_3dDomainInput):
    """3d lid driven cavity, moving wall: z+ wall."""
    def __init__(self, l=1, w=1, h=1, domain_name="Lid-Driven-Cavity"):
        assert l > 0 and w > 0 and h > 0, f"l, w, h = {l}, {w}, {h} is wrong."

        self._lwh_ = [l, w, h]

        super().__init__(domain_name=domain_name)

        x0 = 0
        x1 = l
        y0 = 0
        y1 = w
        z0 = 0
        z1 = h

        self.region_corner_coordinates = {'R:R': ((x0, y0, z0), (x1, y0, z0), (x0, y1, z0), (x1, y1, z0),
                                                  (x0, y0, z1), (x1, y0, z1), (x0, y1, z1), (x1, y1, z1))}
        self.region_side_types = dict() # all plane
        self.boundary_region_sides = {'WALL': ("R:R-N", "R:R-S", "R:R-W", "R:R-E", "R:R-B"),
                                      'LID': "R:R-F"}
        self.region_interpolators = {'R:R': 'transfinite'}
        self.region_type_wr2_metric = {'R:R': 'transfinite'}
        self.internal_parameters = list()  # has to be defined after the super().__init__

    @property
    def lwh(self):
        return self._lwh_