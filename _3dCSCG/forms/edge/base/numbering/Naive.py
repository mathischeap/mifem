# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('./')

from screws.freeze.main import FrozenOnly
from tools.linear_algebra.gathering.chain_matrix.main import Gathering_Matrix, Gathering_Vector



class _3dCSCG_Edge_Numbering_Naive(FrozenOnly):
    def __init__(self, ef):
        self._ef_ = ef
        self._mesh_ = ef.mesh
        self._freeze_self_()

    def _3dCSCG_0Edge(self):
        """Do the numbering if it is an edge 0-form:
        :class:`_3dCSCG.form.edge._0_edge._0Edge`.

        :returns: A tuple of 4 outputs:

            1. (Gathering_Matrix) -- The global numbering in local elements.
            2. (dict) -- The global numbering in local edge elements.
            3. (int) -- Number of dofs in this core.
            4. (None,...) -- Extra numbering information.
        """
        if self._ef_.numbering._parameters_ == dict():
            return self._0Edge_no_parameters()
        else:
            raise NotImplementedError()

    def _0Edge_no_parameters(self):
        """
        Do the numbering if it is an edge 0-form:
        :class:`_3dCSCG.form.edge._0_edge._0Edge`.

        :returns: A tuple of 4 outputs:

            1. (Gathering_Matrix) -- The global numbering in local elements.
            2. (dict) -- The global numbering in local edge elements. Edge-element-wise.
            3. (int) -- Number of dofs in this core.
            4. (None,...) -- Extra numbering information.
        """
        GM = dict()
        GM_EEW = dict()
        local_num_dofs = 0
        extraInfo = None

        _p_ = self._ef_.space.p
        p = [_p_[0] + 1, _p_[1] + 1, _p_[2] + 1]
        tAn = self._mesh_.edge.elements.___PRIVATE_find_type_and_amount_numbered_before___()
        _D_ = {'NS':p[0], 'WE':p[1], 'BF':p[2]}

        for i in self._mesh_.edge.elements:
            e_e_i = self._mesh_.edge.elements[i]
            am_NS, am_WE, am_BF = tAn[i]
            start_num = am_NS * p[0] + am_WE * p[1] + am_BF * p[2]
            GM_EEW[i] = range(start_num, start_num + _D_[e_e_i.direction])
            local_num_dofs += _D_[e_e_i.direction]

        MAP = self._mesh_.edge.elements.map
        for i in MAP:
            vector = tuple()
            for t_j in MAP[i]:
                vector += (GM_EEW[t_j],)
            GM[i] = Gathering_Vector(i, vector)
        GM = Gathering_Matrix(GM, mesh_type='_3dCSCG')

        for i in GM_EEW:
            GM_EEW[i] = Gathering_Vector(i, GM_EEW[i])

        return GM, GM_EEW, local_num_dofs, extraInfo








    def _3dCSCG_1Edge(self):
        """Do the numbering if it is an edge 1-form:
        :class:`_3dCSCG.form.edge._1_edge._1Edge`.

        :returns: A tuple of 4 outputs:

            1. (Gathering_Matrix) -- The global numbering in local elements.
            2. (dict) -- The global numbering in local edge elements.
            3. (int) -- Number of dofs in this core.
            4. (None,...) -- Extra numbering information.
        """
        if self._ef_.numbering._parameters_ == dict():
            return self._1Edge_no_parameters()
        else:
            raise NotImplementedError()

    def _1Edge_no_parameters(self):
        """
        Do the numbering if it is an edge 1-form:
        :class:`_3dCSCG.form.edge._1_edge._1Edge`.

        :returns: A tuple of 4 outputs:

            1. (Gathering_Matrix) -- The global numbering in local elements.
            2. (dict) -- The global numbering in local edge elements. Edge-element-wise.
            3. (int) -- Number of dofs in this core.
            4. (None,...) -- Extra numbering information.
        """
        GM = dict()
        GM_EEW = dict()
        local_num_dofs = 0
        extraInfo = None

        p = self._ef_.space.p
        tAn = self._mesh_.edge.elements.___PRIVATE_find_type_and_amount_numbered_before___()
        _D_ = {'NS':p[0], 'WE':p[1], 'BF':p[2]}

        for i in self._mesh_.edge.elements:
            e_e_i = self._mesh_.edge.elements[i]
            am_NS, am_WE, am_BF = tAn[i]
            start_num = am_NS * p[0] + am_WE * p[1] + am_BF * p[2]
            GM_EEW[i] = range(start_num, start_num + _D_[e_e_i.direction])
            local_num_dofs += _D_[e_e_i.direction]

        MAP = self._mesh_.edge.elements.map
        for i in MAP:
            vector = tuple()
            for t_j in MAP[i]:
                vector += (GM_EEW[t_j],)
            GM[i] = Gathering_Vector(i, vector)
        GM = Gathering_Matrix(GM, mesh_type='_3dCSCG')

        for i in GM_EEW:
            GM_EEW[i] = Gathering_Vector(i, GM_EEW[i])

        return GM, GM_EEW, local_num_dofs, extraInfo









if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\form\edge\numbering\Naive.py
    from _3dCSCG.master import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('crazy', c=0.25)([3,4,5])
    space = SpaceInvoker('polynomials')([('Lobatto',2), ('Lobatto',1), ('Lobatto',3)])
    FC = FormCaller(mesh, space)

    e = FC('0-e')

    print(e.numbering.gathering)