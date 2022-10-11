# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/30 4:58 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.miUsGrid.triangular.forms.standard._1.base.IDT import miUs_Triangular_S1F_InterfaceDofTopology


class miUs_Triangular_oS1F_InterfaceDofTopology(miUs_Triangular_S1F_InterfaceDofTopology):
    """"""

    def __init__(self, sf):
        """"""
        super(miUs_Triangular_oS1F_InterfaceDofTopology, self).__init__(sf)

    @property
    def EEs(self):
        """The `Element-Edges` whose dofs (as well as their basis functions) will receive a minus."""
        mesh = self._sf_.mesh
        eMap = mesh.elements.map
        _EEs_ = dict()
        for e in eMap:
            edge_maps = eMap[e]
            for j, emp in enumerate(edge_maps):
                if emp.isalnum():  # a boundary element edge
                    pass
                else:
                    other = int(emp[:-2])
                    if other < e:  # we put all 'minus' to the element of larger number.
                        other_edge = emp[-1]
                        self_edge = str(j)

                        if other_edge == self_edge:
                            if e not in _EEs_:
                                _EEs_[e] = list()
                            _EEs_[e].append(j)

                        else:
                            if j != 1 and other_edge != '1':
                                if e not in _EEs_:
                                    _EEs_[e] = list()
                                _EEs_[e].append(j)

        return _EEs_



if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
