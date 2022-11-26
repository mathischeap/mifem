# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly


class _2dCSCG_FormsAllocator(FrozenOnly):
    """"""


    @classmethod
    def ___form_name___(cls):
        return {'0-f-i': "_2dCSCG_0Form_Inner", # d on inner 0-form is grad
                '1-f-i': "_2dCSCG_1Form_Inner", # d on inner 1-form is rot
                '2-f-i': "_2dCSCG_2Form_Inner",

                '0-f-o': "_2dCSCG_0Form_Outer", # d on outer 0-form is curl
                '1-f-o': "_2dCSCG_1Form_Outer", # d on outer 1-form is div
                '2-f-o': "_2dCSCG_2Form_Outer",

                '1-t-o': "_2dCSCG_1Trace_Outer",
                }


    @classmethod
    def ___form_path___(cls):
        base_path = '.'.join(str(cls).split(' ')[1][1:-2].split('.')[:-2]) + '.'
        return {'0-f-i': base_path + "standard._0_form.inner.main", # d on inner 0-form is grad
                '1-f-i': base_path + "standard._1_form.inner.main", # d on inner 1-form is rot
                '2-f-i': base_path + "standard._2_form.inner.main",

                '0-f-o': base_path + "standard._0_form.outer.main", # d on outer 0-form is curl
                '1-f-o': base_path + "standard._1_form.outer.main", # d on outer 1-form is div
                '2-f-o': base_path + "standard._2_form.outer.main",

                '1-t-o': base_path + "trace._1_trace.outer.main",
                }