# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/28/2022 12:25 PM
"""
from components.freeze.main import FrozenOnly

from importlib import import_module


class _3dCSCG_DiscreteFieldsAllocator(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        assert mesh.__class__.__name__ == "_3dCSCG_Mesh", f"I need a 3d CSCG mesh."
        self._mesh_ = mesh
        self._freeze_self_()

    def __call__(self, ID, **kwargs):
        """"""
        cls_name = self.___discrete_field_path___()[ID]
        cls_path = self.___discrete_field_name___()[ID]
        clas_body = getattr(import_module(cls_path), cls_name)
        return clas_body(self._mesh_, **kwargs)

    @classmethod
    def ___discrete_field_name___(cls):
        """"""
        return {
            "scalar": "_3dCSCG_DF_Scalar",
            "vector": "_3dCSCG_DF_Vector",
        }

    @classmethod
    def ___discrete_field_path___(cls):
        """"""
        base_path = '.'.join(str(cls).split(' ')[1][1:-2].split('.')[:-2]) + '.'

        return {
            "scalar": base_path + "scalar.main",
            "vector": base_path + "vector.main",
        }
