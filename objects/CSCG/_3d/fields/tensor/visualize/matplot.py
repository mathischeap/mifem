# -*- coding: utf-8 -*-
from screws.freeze.main import FrozenOnly

# import matplotlib.pyplot as plt
# from matplotlib import cm
# from root.config.main import *

class _3dCSCG_TensorField_matplot_Visualize(FrozenOnly):
    """"""
    def __init__(self, f):
        """"""
        self._f_ = f
        self._freeze_self_()


    def __call__(self, *args, **kwargs):
        """"""
        return self.boundary_values(*args, **kwargs)

    def boundary_values(self, *args, **kwargs):
        """"""
        raise NotImplementedError(f"can not matplot boundary values of {self._f_}.")
