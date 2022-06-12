# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/17 7:52 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.mpRfT.base.mesh.main import mpRfT_MeshBase
from objects.mpRfT._2d.mesh.basic_cells.main import mpRfT2_Mesh_BasicCells
from objects.mpRfT._2d.mesh.do.main import mpRfT2_Mesh_Do
from objects.mpRfT._2d.mesh.rcWds.main import mpRfT2_Mesh_RootCellWiseDataStructure
from objects.mpRfT._2d.mesh.refinements.main import mpRfT2_Mesh_Refinements
from objects.mpRfT._2d.mesh.visualize import mpRfT2_Mesh_Visualize
from objects.mpRfT._2d.mesh.coo_map.main import mpRfT2_Mesh_CooMap
from objects.mpRfT._2d.mesh.space.main import mpRfT2_Mesh_Space
from objects.mpRfT._2d.mesh.segments.main import mpRfT2_Mesh_Segments
from objects.mpRfT._2d.mesh.all_root_cells import mpRfT2_Mesh_AllRootCells
from objects.mpRfT._2d.mesh.rcMetricComputing.main import mpRfT2_Mesh_rcMC
from objects.mpRfT._2d.mesh.sgMetricComputing.main import mpRfT2_Mesh_sgMC
from objects.mpRfT._2d.mesh.boundaries.main import mpRfT2_Mesh_Boundaries



class mpRfT2_Mesh(mpRfT_MeshBase):
    """"""

    def __init__(self, cscg, dN, rfd, SNM=1):
        """

        Parameters
        ----------
        cscg
        dN
        rfd : dict
            Refinements dict.
        SNM :
            Segment degree (N) Method.

        """
        super(mpRfT2_Mesh, self).__init__(cscg, dN, rfd)
        self._do_ = mpRfT2_Mesh_Do(self)
        self._rcWds_ = None
        self._refinements_ = None
        self._visualize_ = None
        self._coo_map_ = None
        self._space_ = None
        self._ARC_ = None
        self._rcMC_ = mpRfT2_Mesh_rcMC(self)
        self._sgMC_ = mpRfT2_Mesh_sgMC(self)
        self._segments_ = mpRfT2_Mesh_Segments(self)
        self._segments_.___Pr_find_N_for_all_segments___(SNM)
        self._boundaries_ = mpRfT2_Mesh_Boundaries(self)
        self._freeze_self_()



    def ___Pr_make_basic_cells___(self, dN, rfd, _):
        """

        Parameters
        ----------
        dN
        rfd : dict
            {__repr__ : N, ...} : keys repr of a root-cell, N the space degree of the root-cell.

            Therefore, all sub-cells must be set `N`.
        _

        Returns
        -------

        """
        basic_cells = mpRfT2_Mesh_BasicCells(self) # change this for 3-d mpRfT3_Mesh
        return super(mpRfT2_Mesh, self).___Pr_make_basic_cells___(dN, rfd, basic_cells)

    @property
    def do(self):
        return self._do_

    @property
    def rcWds(self):
        if self._rcWds_ is None:
            self._rcWds_ = mpRfT2_Mesh_RootCellWiseDataStructure(self)
        return self._rcWds_

    @property
    def refinements(self):
        if self._refinements_ is None:
            self._refinements_ = mpRfT2_Mesh_Refinements(self)
        return self._refinements_

    @property
    def visualization(self):
        if self._visualize_ is None:
            self._visualize_ = mpRfT2_Mesh_Visualize(self)
        return self._visualize_

    @property
    def coo_map(self):
        """"""
        if self._coo_map_ is None:
            self._coo_map_ = mpRfT2_Mesh_CooMap(self)
        return self._coo_map_

    @property
    def space(self):
        """"""
        if self._space_ is None:
            self._space_ = mpRfT2_Mesh_Space(self)
        return self._space_

    @property
    def segments(self):
        """"""
        return self._segments_

    @property
    def arc(self):
        """"""
        if self._ARC_ is None:
            self._ARC_ = mpRfT2_Mesh_AllRootCells(self)
        return self._ARC_

    @property
    def rcMC(self):
        """"""
        return self._rcMC_

    @property
    def sgMC(self):
        """"""
        return self._sgMC_

    @property
    def boundaries(self):
        """"""
        return self._boundaries_





if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/mesh/main.py

    # from objects.mpRfT._2d.master import MeshGenerator
    # mesh = MeshGenerator('rectangle')([3,3], 2, show_info=True)

    from __init__ import rfT2
    mesh = rfT2.rm(10, refinement_intensity=0.5)
    # for rp in mesh.arc:
    #     mesh.visualize(rp=rp)

    mesh.visualization()