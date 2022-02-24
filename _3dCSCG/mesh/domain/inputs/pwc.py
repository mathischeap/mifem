



from _3dCSCG.mesh.domain.inputs.base import _3dDomainInputBase




class Parallel_Wall_Channel(_3dDomainInputBase):
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

