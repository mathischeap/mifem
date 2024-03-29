# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/10/13 9:19 PM
"""
from components.freeze.main import FrozenOnly
from objects.CSCG.base.forms.base.BC.interpret.local import CSCG_FORM_BC_Interpret_Local
from objects.CSCG.base.forms.base.BC.interpret.helpers.boundary_cochain import CSCG_FORM_BC_Interpret_BoundaryCochain
from objects.CSCG.base.forms.base.BC.interpret.helpers.boundary_integration import \
    CSCG_FORM_BC_Interpret_BoundaryIntegration
from tools.elementwiseCache.dataStructures.objects.columnVector.main import EWC_ColumnVector


class CSCG_FORM_BC_Interpret(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._mesh_ = f.mesh
        self._bc_ = None
        self._freeze_self_()

    @property
    def local(self):
        """We interpret the BC in local elements."""
        return CSCG_FORM_BC_Interpret_Local(self._f_)
    
    @property
    def boundary_cochain(self):
        """"""
        # Do not cache it. It will change in realtime anyway according to `BC.CF` and `BC.boundaries`
        # because we may apply outer customization to it.
        return EWC_ColumnVector(
            self._f_.mesh.elements,
            CSCG_FORM_BC_Interpret_BoundaryCochain(self._f_),
            'no_cache',
        )

    @property
    def boundary_integration(self):
        # Do not cache it. It will change in realtime anyway according to `BC.CF` and `BC.boundaries`
        # because we may apply outer customization to it.
        return EWC_ColumnVector(
            self._f_.mesh.elements,
            CSCG_FORM_BC_Interpret_BoundaryIntegration(self._f_),
            'no_cache',
        )
