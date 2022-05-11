# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/07 3:06 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.nCSCG.rf2._2d.mesh.cell.frame.segments.main import FrameSegments




class _2nCSCG_CellFrame(FrozenOnly):
    """"""
    def __init__(self, cell):
        """"""
        assert cell.mesh._locker_, f"mesh should be locked before accessing frame. Update mesh first."
        self._cell_ = cell
        self._U = None
        self._D = None
        self._L = None
        self._R = None
        self._freeze_self_()

    def ___Pr_make_Segments___(self):
        """"""
        cell = self._cell_
        i = cell.indices[0]

        ASG = cell.mesh.segments.BCW[i]

        o_and_d = cell.coordinate_transformation.origin_and_delta
        o, delta = o_and_d
        ox, oy = o
        ex, ey = ox + delta, oy + delta

        ySignature = ''
        for _ in cell.indices[1:]:
            if _ in (0, 1):
                ySignature += '0'
            elif _ in (2, 3):
                ySignature += '1'
            else:
                raise Exception()

        xSignature = ''
        for _ in cell.indices[1:]:
            if _ in (0, 2):
                xSignature += '0'
            elif _ in (1, 3):
                xSignature += '1'
            else:
                raise Exception()

        atb = cell.IS.attached_to_Lv0cell_boundary

        USG = dict()
        DSG = dict()
        LSG = dict()
        RSG = dict()

        SG_ci = 'SG-c' + str(i) + ':'
        Ui_signature = SG_ci + 'x' + str(ox) + ':y' + ySignature
        Di_signature = SG_ci + 'x' + str(ex) + ':y' + ySignature
        Li_signature = SG_ci + 'y' + str(oy) + ':x' + xSignature
        Ri_signature = SG_ci + 'y' + str(ey) + ':x' + xSignature

        if atb:
            atUb = cell.IS.attached_to_Lv0cell_U_boundary
            atDb = cell.IS.attached_to_Lv0cell_D_boundary
            atLb = cell.IS.attached_to_Lv0cell_L_boundary
            atRb = cell.IS.attached_to_Lv0cell_R_boundary

            Tmap = cell.mesh.Lv0trace.elements.map[i]
            bt0, bt1, bt2, bt3 = ['t'+str(_) for _ in Tmap]
            Ub_signature = 'SG-' + bt0 + '-y' + ySignature
            Db_signature = 'SG-' + bt1 + '-y' + ySignature
            Lb_signature = 'SG-' + bt2 + '-x' + xSignature
            Rb_signature = 'SG-' + bt3 + '-x' + xSignature

            for _r in ASG:
                segment = ASG[_r]

                nfiYet = True # not find it yet!

                #---- looking for U segments -------
                if atUb: # on the base-cell Upper boundary
                    if Ub_signature in _r: # we find a lv0-trace-segment
                        USG[_r] = segment
                        nfiYet = False
                else:
                    if Ui_signature in _r:
                        USG[_r] = segment
                        nfiYet = False

                if nfiYet:
                    #---- looking for D segments -------
                    if atDb: # on the base-cell Down boundary
                        if Db_signature in _r: # we find a lv0-trace-segment
                            DSG[_r] = segment
                            nfiYet = False
                    else:
                        if Di_signature in _r:
                            DSG[_r] = segment
                            nfiYet = False

                    if nfiYet:
                        #---- looking for L segments -------
                        if atLb: # on the base-cell Left boundary
                            if Lb_signature in _r: # we find a lv0-trace-segment
                                LSG[_r] = segment
                                nfiYet = False
                        else:
                            if Li_signature in _r:
                                LSG[_r] = segment
                                nfiYet = False

                        if nfiYet:
                            #---- looking for R segments -------
                            if atRb: # on the base-cell Left boundary
                                if Rb_signature in _r: # we find a lv0-trace-segment
                                    RSG[_r] = segment
                            else:
                                if Ri_signature in _r:
                                    RSG[_r] = segment

        else:
            for _r in ASG:
                if Ui_signature in _r:
                    USG[_r] = ASG[_r]
                elif Di_signature in _r:
                    DSG[_r] = ASG[_r]
                elif Li_signature in _r:
                    LSG[_r] = ASG[_r]
                elif Ri_signature in _r:
                    RSG[_r] = ASG[_r]
                else:
                    pass

        self._U = FrameSegments(cell, 'U', USG)
        self._D = FrameSegments(cell, 'D', DSG)
        self._L = FrameSegments(cell, 'L', LSG)
        self._R = FrameSegments(cell, 'R', RSG)

    @property
    def U(self):
        if self._U is None:
            self.___Pr_make_Segments___()
        return  self._U

    @property
    def D(self):
        if self._D is None:
            self.___Pr_make_Segments___()
        return  self._D

    @property
    def L(self):
        if self._L is None:
            self.___Pr_make_Segments___()
        return  self._L

    @property
    def R(self):
        if self._R is None:
            self.___Pr_make_Segments___()
        return  self._R






if __name__ == "__main__":
    # mpiexec -n 4 python objects/nCSCG/rf2/_2d/mesh/cell/frame/main.py
    from objects.nCSCG.rf2._2d.__tests__.Random.mesh import random_mesh_of_elements_around as rm2
    mesh = rm2(100, refinement_intensity=0.5)

    # from objects.nCSCG.rf2._2d.master import MeshGenerator
    # mesh = MeshGenerator('crazy')([2, 3], EDM='chaotic', show_info=True)
    # mesh.do.unlock()
    #
    # if 0 in mesh.cscg.elements:
    #     c0 = mesh(0)
    #     c0.do.refine()
    #     c00 = c0(0)
    #     c00.do.refine()
    #     C = c0(2)
    #
    # mesh.do.update()
    #
    # if 0 in mesh.cscg.elements:
    #     D = C.frame.D
    #     for segment in D:
    #         print(segment.where)

    for ind in mesh:
        F = mesh(ind).frame
        D= F.D