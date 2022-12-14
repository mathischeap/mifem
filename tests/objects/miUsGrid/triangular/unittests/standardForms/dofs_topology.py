# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 9/29/2022 10:19 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
import numpy as np

from components.freeze.base import FrozenOnly
from components.miscellaneous.miprint import miprint
from __init__ import miTri

from root.config.main import RANK, MASTER_RANK, COMM

def u_fun(t, x, y): return np.pi * np.exp(np.pi * x) * np.sin(np.pi * y) + 0 * t
def v_fun(t, x, y): return np.pi * np.sin(np.pi * x) * np.cos(0.983*np.pi * y) + 0 * t

class Test_dofs_topology_S1F(FrozenOnly):
    """"""

    def __init__(self):
        """"""
        miprint("miTopo1Dofs -> [Test_dofs_topology_S1F] ...")
        self.fc = miTri.call('rand0', 2)
        self.mesh = self.fc.mesh
        self._freeze_self_()

    def __call__(self):

        vector = self.fc('vector', [u_fun, v_fun])
        vector.current_time = 0

        eMap = self.mesh.elements.map

        for f in [self.fc('1-f-o'), ]:

            f.CF = vector
            f.discretize()
            gm = f.numbering.gathering

            dof_cochain = dict()
            for e in self.mesh.elements:
                local_cochain = f.cochain.local[e]
                full_vector = gm[e].full_vector
                for local_number, global_number in enumerate(full_vector):
                    if global_number not in dof_cochain:
                        dof_cochain[global_number] = [local_cochain[local_number],]
                    else:
                        dof_cochain[global_number].append(local_cochain[local_number])

            dof_cochain = COMM.gather(dof_cochain, root=MASTER_RANK)

            if RANK == MASTER_RANK:
                DOF_COCHAIN = dict()
                for dc in dof_cochain:
                    for dof in dc:
                        if dof in DOF_COCHAIN:
                            DOF_COCHAIN[dof].extend(dc[dof])
                        else:
                            DOF_COCHAIN[dof] = dc[dof]

                for dof in DOF_COCHAIN:
                    if len(DOF_COCHAIN[dof]) == 2:
                        d0, d1 = DOF_COCHAIN[dof]
                        if abs(d0 - d1) < 1e-10:
                            pass
                        else:
                            assert abs(d0 + d1) < 1e-10
                    else:
                        assert len(DOF_COCHAIN[dof]) == 1

            dof_map_dict = dict()
            for e in eMap:
                for edge in range(3):
                    dofs = f.numbering.do.find.dofs_on_element_edge(e, edge)
                    for dof in dofs:
                        if eMap[e][edge].isalnum():
                            pass
                        else:
                            if dof not in dof_map_dict:
                                dof_map_dict[dof] = [eMap[e][edge],]
                            else:
                                dof_map_dict[dof].append(eMap[e][edge])

            dof_map_dict = COMM.gather(dof_map_dict, root=MASTER_RANK)
            if RANK == MASTER_RANK:
                DICT = dict()
                for dmd in dof_map_dict:
                    for dof in dmd:
                        if dof not in DICT:
                            DICT[dof] = dmd[dof]
                        else:
                            DICT[dof].extend(dmd[dof])

                for dof in DICT:
                    # noinspection PyUnboundLocalVariable
                    assert len(DICT[dof]) == len(DOF_COCHAIN[dof]) == 2
                    map0, map1 = DICT[dof]
                    d0, d1 = DOF_COCHAIN[dof]
                    if '+' in map0:
                        assert '+' in map1
                    else:
                        assert '-' in map0 and '-' in map1

                    assert abs(d0 - d1) < 1e-10, f"{DICT[dof], DOF_COCHAIN[dof]}"

        scalar = self.fc('scalar', u_fun)
        scalar.current_time = 0
        for f in [self.fc('0-f-o'), self.fc('0-f-i')]:

            f.CF = scalar
            f.discretize()
            gm = f.numbering.gathering

            dof_cochain = dict()
            for e in self.mesh.elements:
                local_cochain = f.cochain.local[e]
                full_vector = gm[e].full_vector
                for local_number, global_number in enumerate(full_vector):
                    if global_number not in dof_cochain:
                        dof_cochain[global_number] = [local_cochain[local_number],]
                    else:
                        dof_cochain[global_number].append(local_cochain[local_number])

            dof_cochain = COMM.gather(dof_cochain, root=MASTER_RANK)

            if RANK == MASTER_RANK:
                DOF_COCHAIN = dict()
                for dc in dof_cochain:
                    for dof in dc:
                        if dof in DOF_COCHAIN:
                            DOF_COCHAIN[dof].extend(dc[dof])
                        else:
                            DOF_COCHAIN[dof] = dc[dof]

                for dof in DOF_COCHAIN:
                    if len(DOF_COCHAIN[dof]) > 1:
                        d0 = DOF_COCHAIN[dof][0]
                        for d in DOF_COCHAIN[dof][1:]:
                            assert abs(d - d0) < 1e-10
                    else:
                        assert len(DOF_COCHAIN[dof]) == 1

        return 1

if __name__ == '__main__':
    # mpiexec -n 4 python tests/objects/miUsGrid/triangular/unittests/standardForms/dofs_topology.py
    Test_dofs_topology_S1F()()
