# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11:10 PM
"""
import sys
if './' not in sys.path: sys.path.append('./')

from objects.nCSCG.rf2._2d.mesh.space.base.base import _2nCSCG_SpaceBase
from objects.nCSCG.rf2._2d.mesh.space.polynomials.do.main import _2nCSCG_SpacePolyDo
from objects.CSCG._2d.spaces.polynomials import _2dCSCG_PolynomialSpace

import numpy as np



class _2nCSCG_PolynomialSpace(_2nCSCG_SpaceBase):
    """"""
    def __init__(self, dN, nodes_type='Lobatto'):
        """

        Parameters
        ----------
        dN : int
        nodes_type : str
            {'Lobatto', }

        """
        super(_2nCSCG_PolynomialSpace, self).__init__(dN)
        self._do_ = None
        self._visualize_ = None
        self.___node_type___ = nodes_type
        self._r = None
        self._GoN_ = dict() # gird_of_nodes
        self._GoN_ravel_ = dict() # gird_of_nodes raveled.

        self._freeze_self_()
        #----- customize following parameters for particular spaces -------------------------------------
        self._PRM = (dN, 'polynomials', {'nodes_type':nodes_type})
        self.___space_class___ = _2dCSCG_PolynomialSpace
        #================================================================================================
        self.___pool___[self.dN] = self.___space_class___((self.___node_type___, dN), ndim=2)
        # so we only use Nx = Ny.

        dNodes = self.___pool___[self.dN].nodes
        GoN = np.meshgrid(*dNodes, indexing='ij')
        self._GoN_[self.dN] = GoN
        # noinspection PyUnresolvedReferences
        self._GoN_ravel_[self.dN] = GoN[0].ravel('F'), GoN[1].ravel('F')

    def __repr__(self):
        """"""
        if self._r is None:
            r = self.__class__.__name__ + ':' + self.___node_type___ + '|' + str(self._dN_) + '-'
            Ns = list(self.___pool___.keys())
            Ns.sort()
            r += str(Ns)
            self._r = r
        return self._r

    @property
    def do(self):
        """"""
        if self._do_ is None:
            self._do_ = _2nCSCG_SpacePolyDo(self)
        return self._do_

    @property
    def GoN(self):
        return self._GoN_

    @property
    def GoN_ravel(self):
        return self._GoN_ravel_




if __name__ == '__main__':
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/mesh/space/polynomials/main.py
    from objects.nCSCG.rf2._2d.__tests__.Random.mesh import random_mesh_of_elements_around as rm2
    mesh = rm2(100, refinement_intensity=0.5)
    space = _2nCSCG_PolynomialSpace(3)