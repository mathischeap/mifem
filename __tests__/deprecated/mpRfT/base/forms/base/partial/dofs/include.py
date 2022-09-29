# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/13 9:11 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class mpRfT_Form_PartialDofs_Include(FrozenOnly):
    """"""

    def __init__(self, pd):
        """"""
        self._f_ = pd._f_
        self._mesh_ = pd._mesh_
        self._pd_ = pd

        self._freeze_self_()

    def boundaries(self, boundary_names):
        """We will add keys to self._pd_.keys. Those added keys represents the surfaces of root-cells
        on the boundary names.
        """
        if isinstance(boundary_names, str):
            boundary_names = [boundary_names,]
        bns = self._f_.mesh.boundaries.names
        for bn in boundary_names:
            assert bn in bns, f"valid boundary {bn} is not a boundary name."

        rrc = self._mesh_.boundaries.rrc
        indicators = self._pd_._indicators_
        for bn in boundary_names:
            root_cells_and_edges = rrc[bn]
            for rce in root_cells_and_edges:
                rc, edge = rce
                if rc in indicators:
                    indicators[rc].append(edge)
                else:
                    indicators[rc] = [edge,]






if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/base/forms/base/partial/dofs/include.py
    from __init__ import rfT2

    fc = rfT2.rf(100, mesh_pool='crazy')

    f1 = fc('1-f-o')
    f1.BC.valid_boundaries = ['Upper', 'Left']
    fpd = f1.BC.partial_dofs

    t1 = fc('nst')
    t1.BC.valid_boundaries = ['Down', 'Right']
    tpd = t1.BC.partial_dofs


    print(tpd._indicators_)