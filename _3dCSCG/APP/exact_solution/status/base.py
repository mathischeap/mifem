# -*- coding: utf-8 -*-


from root.config.main import *
from screws.quadrature import Quadrature
from screws.freeze.main import FrozenOnly



class Base(FrozenOnly):
    """
    A base (parent) for all exact solution classes.
    """
    def __init__(self, es):
        self._es_ = es
        self._mesh_ = es.mesh

    # ... valid time: override it to make the exact solution only valid at certain time. ...
    @property
    def valid_time(self):
        return None

    @property
    def mesh(self):
        """
        all exact solutions have mesh, no matter if mesh is useful (general exact solutions not use mesh).

        :return: A ``_3dCSCG_Mesh`` object.
        """
        return self._mesh_