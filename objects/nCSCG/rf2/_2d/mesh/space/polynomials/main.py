# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11:10 PM
"""
import sys
if './' not in sys.path: sys.path.append('./')

from objects.nCSCG.rf2._2d.mesh.space.base import _2nCSCG_SpaceBase
from objects.nCSCG.rf2._2d.mesh.space.polynomials.do import _2nCSCG_SpacePolyDo
from objects.nCSCG.rf2._2d.mesh.space.polynomials.visualize import _2nCSCG_PolySpaceVisualize
from objects.CSCG._2d.spaces.polynomials import _2dCSCG_PolynomialSpace




class _2nCSCG_PolynomialSpace(_2nCSCG_SpaceBase):
    """"""
    def __init__(self, N, nodes_type='Lobatto'):
        """

        Parameters
        ----------
        N : int
        nodes_type : str
            {'Lobatto', }

        """
        super(_2nCSCG_PolynomialSpace, self).__init__(N)
        self._do_ = None
        self._visualize_ = None
        self.___node_type___ = nodes_type
        self._r = None
        self._freeze_self_()
        #----- customize following parameters for particular spaces -------------------------------------
        self._PRM = (N, 'polynomials', {'nodes_type':nodes_type})
        self.___space_class___ = _2dCSCG_PolynomialSpace
        #================================================================================================
        self.___pool___[self._dN_] = self.___space_class___((self.___node_type___, N), ndim=2) # so we only use Nx = Ny.

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
    def visualize(self):
        """"""
        if self._visualize_ is None:
            self._visualize_ = _2nCSCG_PolySpaceVisualize(self)
        return self._visualize_




if __name__ == '__main__':
    # mpiexec -n 4 python objects/nCSCG/rf2/_2d/mesh/space/polynomials/main.py
    from objects.nCSCG.rf2._2d.__tests__.Random.mesh import random_mesh_of_elements_around as rm2
    mesh = rm2(100, refinement_intensity=0.5)
    space = _2nCSCG_PolynomialSpace(3)