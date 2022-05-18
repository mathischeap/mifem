# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/17 2:51 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.nCSCG.rf2._2d.mesh.refinement.p_.base import _2nCSCG_RF2_p_Refinement_Base


class pR_Simple(_2nCSCG_RF2_p_Refinement_Base):
    """A simple p-refinement is one that only for the mesh with no h-refinement at all."""

    def __init__(self, mesh):
        """"""
        super(pR_Simple, self).__init__(mesh)
        assert self.mesh.IS.basic, f"simple p-refinement only works for basic mesh."
        self._freeze_self_()

    # def




if __name__ == "__main__":
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/mesh/refinement/p_/simple.py
    from objects.nCSCG.rf2._2d.master import MeshGenerator

    mesh = MeshGenerator('crazy')([3, 3], 2, EDM='chaotic', show_info=True)
    # mesh.visualize()
    pr = mesh.refinement.p('simple')


