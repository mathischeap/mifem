# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 12/9/2022 3:20 PM
"""
from components.freeze.main import FrozenOnly
from root.config.main import COMM, MPI


class _3dCSCG_TraceForm_DofWhether(FrozenOnly):
    """"""

    def __init__(self, dof):
        """"""
        self._dof_ = dof
        self._on_periodic_boundary_ = None
        self._freeze_self_()

    @property
    def on_periodic_boundary(self):
        """"""
        if self._on_periodic_boundary_ is None:
            te = self._dof_.trace_element_position
            if te is not None:
                te = te[0]
                ToF = self._dof_._tf_.mesh.trace.elements[te].whether.on_periodic_boundary

            else:
                ToF = False

            self._on_periodic_boundary_ = COMM.allreduce(ToF, op=MPI.LOR)

        return self._on_periodic_boundary_
