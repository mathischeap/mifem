# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/22 7:47 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from objects.mpRfT._2d.mesh.coo_map.helpers.base import mpRfT2_CooMapBase
from components.quadrature import Quadrature
import numpy as np



class mpRfT2_Mesh_PartialSegmentIntegral(mpRfT2_CooMapBase):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._segCache_ = dict()
        self._nst_ = None
        self._freeze_self_()

    @property
    def distribution(self):
        return 'PSI'

    @property
    def ___Pr_rcMC_key___(self):
        raise NotImplementedError()

    def ___Pr_rcMC_nodes___(self, rp):
        raise NotImplementedError()

    @property
    def ___Pr_sgMC_key___(self):
        """A key implying the value for metric involved computing in each root-cell."""
        raise NotImplementedError()

    def ___Pr_sgMC_nodes___(self, rp):
        raise NotImplementedError()

    def __call__(self, nst, qdp=0):
        """

        Parameters
        ----------
        nst : mpRfT2_NSgF
        qdp : int

        Returns
        -------
        edge_wise_NW : dict
            {'U': (nodes, -1, weights),
             'D': (nodes,  1, weights),
             'L': (-1, nodes, weights),
             'R': ( 1, nodes, weights)}

            where nodes is of shape (N*2,) and weights are of shape (N,)

        """
        assert nst.__class__.__name__ == 'mpRfT2_NSgF'
        self._nst_ = nst
        return self

    def __getitem__(self, rc_rp):
        """"""
        cell = self._mesh_[rc_rp]
        frame = cell.frame

        edge_wise_NW = dict()

        for edge in frame:

            segments = frame[edge]

            mark = segments.type_wrt_metric

            if isinstance(mark, str):
                mark += edge

            OO, DD = segments.origin_and_delta

            if isinstance(OO, tuple):
                OOx, OOy = OO
                if segments.direction == 'UD':
                    OO = OOx
                else:
                    OO = OOy
            else:
                pass

            if isinstance(mark, str) and mark in self._segCache_:
                e_nw = self._segCache_[mark]

            else:

                e_nw = list()

                for seg in segments:

                    o, d = seg.origin_and_delta
                    if isinstance(o, tuple):
                        ox, oy = o
                        if seg.direction == 'UD':
                            o = ox
                        else:
                            o = oy
                    else:
                        pass

                    N = self._nst_.N[seg]

                    nodes, weights = Quadrature(N).quad
                    intervals = Quadrature(N+1, category='Lobatto').quad[0]
                    diff = np.diff(intervals)

                    NODES = list()
                    for i, start in enumerate(intervals[:-1]):

                        nds = start + (nodes + 1) * diff[i] / 2

                        NODES .append( nds )

                    NODES = np.concatenate(NODES)

                    o = 2 * (o - OO) / DD - 1
                    d = 2 * d / DD
                    NODES = o + (NODES + 1 ) * d / 2

                    if seg.direction == 'UD':
                        xi = NODES
                        if edge == 'L':
                            eta = np.array([-1,])
                        elif edge == 'R':
                            eta = np.array([1,])
                        else:
                            raise Exception()

                    elif seg.direction == 'LR':
                        eta = NODES

                        if edge == 'U':
                            xi = np.array([-1,])
                        elif edge == 'D':
                            xi = np.array([1,])
                        else:
                            raise Exception()

                    else:
                        raise Exception()

                    e_nw.append((xi, eta, weights, diff))

                if isinstance(mark, str):
                    self._segCache_[mark] = e_nw

            edge_wise_NW[edge] = e_nw

        return edge_wise_NW



if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/mesh/coo_map/helpers/partial_segment_integral.py

    from __init__ import rfT2
    mesh = rfT2.rm(50, refinement_intensity=0.5,N_range=(2,4))

    FC = rfT2.form(mesh)
    t = FC('nst', ndp=-1, name='trace')

    coo_map = mesh.coo_map.partial_segment_integral(t)

    for rc_rp in mesh.rcfc:
        a = coo_map[rc_rp]