# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/03 10:05 PM
"""
from objects.CSCG._2d.mesh.domain.inputs.base import DomainInputBase
from root.config.main import COMM, RANK, MASTER_RANK


class TriangleTest(DomainInputBase):
    """"""
    def __init__(self):
        """
        """
        super().__init__(domain_name='TriangleTest')

        if RANK == MASTER_RANK:
            x0 = 0.053
            y0 = 0.08712

            x1 = 1.10125
            y1 = 0.5246

            x2 = 0
            y2 = 1.0985

            COO = [x0, x1, x2, y0, y1, y2]
        else:
            COO = None

        COO = COMM.bcast(COO, root=MASTER_RANK)
        x0, x1, x2, y0, y1, y2 = COO

        # _____________ standard inputs ________________________________________________
        self.region_corner_coordinates = {'R:R': ((x0, y0), (x1, y1), (x2, y2), (x1, y1))}
        self.region_edge_types = {}
        self.boundary_region_edges = {'Upper': "R:R-U", 'Down': "R:R-D", 'Left': "R:R-L", 'Right': "R:R-R"}
        self.region_interpolators = 'transfinite'

        self.region_type_wr2_metric = {'R:R': 'transfinite'}
        self.region_sequence = ('R:R',)
        # ------------------------------------------------------------------------------
