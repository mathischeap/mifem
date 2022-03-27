
from objects.CSCG._2d.mesh.domain.inputs.base import DomainInputBase
from screws.decorators.classproperty.main import classproperty


class Quadrangle(DomainInputBase):
    """ A quadrangle of four straight edges.
    """
    def __init__(self, p_UL=(0,0), p_DL=(1,0), p_UR=(1,1), p_DR=(2,1)):

        super().__init__(domain_name='Quadrangle')

        self.region_corner_coordinates = {'R:R': (p_UL, p_DL, p_UR, p_DR)}
        self.region_edge_types = {}
        self.boundary_region_edges = {'Upper': "R:R-U", 'Down': "R:R-D", 'Left': "R:R-L", 'Right': "R:R-R"}
        self.region_interpolators = 'transfinite'
        self.region_type_wr2_metric = {'R:R': 'transfinite'}
        self.region_sequence = ('R:R',)


    @classproperty
    def statistic(cls):
        return {'periodic': False,
                'region num': 1,
                'mesh boundary num': 4, # the amount of mesh boundaries (instead of domain boundaries)
                }

    @classproperty
    def random_parameters(cls):
        return dict()