# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 5/24/2022 9:02 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from objects.mpRfT._2d.mesh.MetricComputing_rc.helpers.base import mpRfT2_Mesh_rcMC_Base


class mpRfT2_Mesh_rcMC_iJM(mpRfT2_Mesh_rcMC_Base):
    """"""

    def __getitem__(self, rp):
        """"""
        cell = self._mesh_[rp]
        mark = self._key_(rp)
        if mark.isnumeric(): # chaotic cell, compute it in real-time
            nodes = self._nodes_(rp)

            #-----------------------------------------------------------------
            return cell.coordinate_transformation.inverse_Jacobian_matrix(*nodes)
            #===================================================================

        else:
            if mark in self._cache_:
                return self._cache_[mark]
            else:
                nodes = self._nodes_(rp)

                #--------------------------------------------------------------
                RETURN = cell.coordinate_transformation.inverse_Jacobian_matrix(*nodes)
                self._cache_[mark] = RETURN
                return RETURN
                #===================================================================



if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
