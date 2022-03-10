
from _2dCSCG.mesh.domain.inputs.base import DomainInputBase

class Rectangle(DomainInputBase):
    """ A rectangle computational whose edges are parallel to the axes.

    --------------------------------------> y
    |
    |   p_UL                              p_UR
    |    -----------------------------------
    |    |                                 |
    |    |                                 |
    |    |                                 |
    |    |                                 |     width
    |    |                                 |
    |    |                                 |
    |    -----------------------------------
    |   p_DL                              p_DR
    |                 length
    v
    x

    """
    def __init__(self, p_UL=(0,0), p_DL=(1,0), p_UR=(0,1), p_DR=(1,1), region_layout=(2,2)):
        """

        :param p_UL:
        :param p_DL:
        :param p_UR:
        :param p_DR:
        :param region_layout: How many regions along each direction. For example,
            `region_layout` = (2,3)
            we will have 2 regions along x direction, and 3 regions along y direction. So all
            regions are structured-distributed.

        """
        super().__init__(domain_name='Rectangle')
        #---- check region_layout : must be a tuple or list of two positive integers ---------------
        assert len(region_layout) == 2 and \
               (region_layout[0] > 0 and region_layout[0] % 1==0) and \
               (region_layout[1] > 0 and region_layout[1] % 1==0), \
            f"region_layout = {region_layout} wrong!"
        # check four points form a rectangle parallel to the axes ----------------------------------
        assert p_UL[1] == p_DL[1] and p_UL[0] < p_DL[0], f"Left edge not parallel to x axis."
        assert p_UL[0] == p_UR[0] and p_UL[1] < p_UR[1], f"Upper edge not parallel to y axis."
        assert p_UR[1] == p_DR[1] and p_UR[0] < p_DR[0], f"Left edge not parallel to x axis."
        assert p_DL[0] == p_DR[0] and p_DL[1] < p_DR[1], f"Upper edge not parallel to y axis."
        width = p_DL[0] - p_UL[0]
        assert width == p_DR[0] - p_UR[0], f"not a rectangle"
        length =  p_UR[1] - p_UL[1]
        assert length == p_DR[1] - p_DL[1], f"not a rectangle"
        # new we parse the four corners of all regions ---------------------------------------------
        L0, L1 = region_layout
        x_step = width / L0
        y_step = length / L1
        x_start, y_start = p_UL
        Rs = [['' for _ in range(L1)] for _ in range(L0)]
        region_corner_coordinates = dict()
        region_sequence = tuple()
        for i in range(L0):
            for j in range(L1):
                region_name = 'R:' + (i+1)*'i' + '_' + (j+1)*'j'
                x = x_start + i * x_step
                y = y_start + j * y_step
                region_corner_coordinates[region_name] = ((x, y),
                                                          (x + x_step, y),
                                                          (x, y + y_step),
                                                          (x + x_step, y + y_step))
                region_sequence += (region_name,)
                Rs[i][j] = region_name

        boundary_region_edges = dict()
        boundary_region_edges['Upper'] = tuple()
        boundary_region_edges['Down'] = tuple()
        boundary_region_edges['Left'] = tuple()
        boundary_region_edges['Right'] = tuple()
        for j in range(L1):
            boundary_region_edges['Upper'] += (Rs[0][j] + '-U',)
            boundary_region_edges['Down'] += (Rs[-1][j] + '-D',)
        for i in range(L0):
            boundary_region_edges['Left'] += (Rs[i][0] + '-L',)
            boundary_region_edges['Right'] += (Rs[i][-1] + '-R',)

        self.region_corner_coordinates = region_corner_coordinates
        self.region_edge_types = dict() # all straight lines
        self.boundary_region_edges = boundary_region_edges
        self.region_interpolators = 'transfinite'
        self.region_type_wr2_metric = 'transfinite'
        self.region_sequence = region_sequence
