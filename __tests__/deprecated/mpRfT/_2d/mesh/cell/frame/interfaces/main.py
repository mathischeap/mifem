# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/22 12:37 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class mpFfT2_CellFrame_Interfaces(FrozenOnly):
    """The four interfaces on U->D->L->R edges."""

    def __init__(self, frame):
        """"""
        self._frame_ = frame
        self._U, self._D, self._L, self._R = None, None, None, None
        self._freeze_self_()

    def ___Pr_find_the_interfaces_surrounding___(self):
        """"""
        frame = self._frame_
        mesh = frame._cell_._mesh_
        for edge in 'UDLR':
            segments = frame[edge]
            sgs_rp = segments.__repr__().split('|')
            assert sgs_rp[0][-1] == edge

            if len(segments) > 1:
                # these segments must form an interface
                if_rp = 'IF>' + '|'.join(sgs_rp[1:])
                IF = mesh.interfaces[if_rp]
                assert not IF.IS.symmetric
            else:
                sg_rp = sgs_rp[1]
                if_rp = 'IF>' + sg_rp

                if if_rp in mesh.interfaces:
                    IF = mesh.interfaces[if_rp]
                    if IF.IS.on_mesh_boundary:
                        assert not IF.IS.symmetric
                    else:
                        assert IF.IS.symmetric

                else:
                    # This edge is within a larger interface.
                    TF = False

                    for if_rp in mesh.interfaces:
                        IF = mesh.interfaces[if_rp]

                        if len(IF) == 1:
                            pass
                        else:
                            if sg_rp in IF.__repr__():
                                TF = True
                                break

                    assert TF

            setattr(self, f"_{edge}", IF)

    @property
    def U(self):
        if self._U is None:
            self.___Pr_find_the_interfaces_surrounding___()
        return self._U

    @property
    def D(self):
        if self._D is None:
            self.___Pr_find_the_interfaces_surrounding___()
        return self._D

    @property
    def L(self):
        if self._L is None:
            self.___Pr_find_the_interfaces_surrounding___()
        return self._L

    @property
    def R(self):
        if self._R is None:
            self.___Pr_find_the_interfaces_surrounding___()
        return self._R

    def __iter__(self):
        for edge in 'UDLR':
            yield edge

    def __getitem__(self, edge):
        return getattr(self, f"_{edge}")





if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/mesh/cell/frame/interfaces/main.py

    from __init__ import rfT2
    mesh = rfT2.rm(10)

    for rc_rp in mesh.rcfc:
        cell = mesh[rc_rp]

        interfaces = cell.frame.interfaces

        U = interfaces.U
        print(U)
