



from objects.CSCG._3d.mesh.domain.inputs.base import _3dDomainInputBase
from screws.decorators.classproperty.main import classproperty



class Lid_Driven_Cavity(_3dDomainInputBase):
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



    @classproperty
    def statistic(cls):
        raise NotImplementedError()


    @classproperty
    def random_parameters(cls):
        raise NotImplementedError()
