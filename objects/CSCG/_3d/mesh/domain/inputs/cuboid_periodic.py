

from objects.CSCG._3d.mesh.domain.inputs.base import _3dDomainInputBase


from screws.decorators.classproperty.main import classproperty
import random
from root.config.main import rAnk, mAster_rank, cOmm



class CuboidPeriodic(_3dDomainInputBase):
    """"""
    def __init__(self, p_NWB=(0,0,0), width=1, length=1, height=1, region_layout=(2,2,2)):
        """

        Parameters
        ----------
        p_NWB : tuple
            The coordinate of the North-West-Back corner point of the cuboid.
        width : float
            x-width
        length : float
            y-length
        height : float
            z-height
        region_layout : tuple
            If `region_layout = (a, b, c)`, there are `a` regions along x-direction, `b` regions
            along y-direction and `c` regions along z-direction.
        """
        super().__init__(domain_name='CuboidPeriodic')
        #---- check region_layout : must be a tuple or list of two positive integers ------------
        assert len(region_layout) == 3 and \
               (region_layout[0] > 0 and region_layout[0] % 1==0) and \
               (region_layout[1] > 0 and region_layout[1] % 1==0) and \
               (region_layout[2] > 0 and region_layout[2] % 1==0), \
            f"region_layout = {region_layout} wrong!"
        assert width > 0 and length > 0 and height > 0
        # new we parse the 8 corners of all regions ---------------------------------------------
        L0, L1, L2 = region_layout
        x_step = width / L0
        y_step = length / L1
        z_step = height / L2
        x_start, y_start, z_start = p_NWB
        Rs = [[['' for _ in range(L2)] for _ in range(L1)] for _ in range(L0)]
        region_corner_coordinates = dict()
        region_sequence = tuple()

        for k in range(L2):
            for j in range(L1):
                for i in range(L0):
                    region_name = 'R:' + (i+1)*'i' + '_' + (j+1)*'j' + '_' + (k+1)*'k'
                    x = x_start + i * x_step
                    y = y_start + j * y_step
                    z = z_start + k * z_step
                    region_corner_coordinates[region_name] = (
                          (x         , y         , z         ),
                          (x + x_step, y         , z         ),
                          (x         , y + y_step, z         ),
                          (x + x_step, y + y_step, z         ),
                          (x         , y         , z + z_step),
                          (x + x_step, y         , z + z_step),
                          (x         , y + y_step, z + z_step),
                          (x + x_step, y + y_step, z + z_step)
                    )

                    region_sequence += (region_name,)
                    Rs[i][j][k] = region_name

        boundary_region_edges = dict()
        boundary_region_edges['North'] = tuple()
        boundary_region_edges['South'] = tuple()
        boundary_region_edges['West'] = tuple()
        boundary_region_edges['East'] = tuple()
        boundary_region_edges['Back'] = tuple()
        boundary_region_edges['Front'] = tuple()

        for k in range(L2):
            for j in range(L1):
                boundary_region_edges['North'] += (Rs[ 0][j][k] + '-N',)
                boundary_region_edges['South'] += (Rs[-1][j][k] + '-S',)

        for k in range(L2):
            for i in range(L0):
                boundary_region_edges['West'] += (Rs[i][ 0][k] + '-W',)
                boundary_region_edges['East'] += (Rs[i][-1][k] + '-E',)

        for j in range(L1):
            for i in range(L0):
                boundary_region_edges['Back']  += (Rs[i][j][ 0] + '-B',)
                boundary_region_edges['Front'] += (Rs[i][j][-1] + '-F',)

        self.region_corner_coordinates = region_corner_coordinates
        self.region_side_types = dict()
        self.boundary_region_sides = boundary_region_edges
        self.region_interpolators = 'orthogonal'
        self.periodic_boundary_pairs = {'South=North': 'regular',
                                        'West=East'  : 'regular',
                                        'Back=Front' : 'regular'}
        self.region_type_wr2_metric = 'orthogonal'
        self.region_sequence = region_sequence


    @classproperty
    def statistic(cls):
        return {'periodic': True,
                'region num': 'unknown',
                'mesh boundary num': 0, # the amount of mesh boundaries (instead of domain boundaries)
                }

    @classproperty
    def random_parameters(cls):
        if rAnk == mAster_rank:
            rp = {'p_NWB': (random.uniform(-1,1), random.uniform(-1,1), random.uniform(-1,1)),
                  'width':random.uniform(1,3),
                  'length':random.uniform(1,3),
                  'height':random.uniform(1,3),
                  "region_layout": random.sample((1,1,2), 3)
                }
        else:
            rp = None

        return cOmm.bcast(rp, root=mAster_rank)