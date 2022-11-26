# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 5/24/2022 9:05 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from objects.mpRfT._2d.mesh.MetricComputing_sg.helpers.base import mpRfT2_Mesh_sgMC_Base


class mpRfT2_Mesh_sgMC_Jacobian(mpRfT2_Mesh_sgMC_Base):
    """"""

    def __getitem__(self, sg):
        """"""
        mark = self._key_(sg)
        if mark.isnumeric(): # chaotic cell, compute it in real-time
            nodes = self._nodes_(sg)

            #-----------------------------------------------------------------
            return sg.coordinate_transformation.Jacobian(*nodes)
            #===================================================================

        else:
            if mark in self._cache_:
                return self._cache_[mark]
            else:
                nodes = self._nodes_(sg)

                #--------------------------------------------------------------
                RETURN = sg.coordinate_transformation.Jacobian(*nodes)
                self._cache_[mark] = RETURN
                return RETURN
                #===================================================================



if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
