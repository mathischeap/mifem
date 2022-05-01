# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from screws.freeze.main import FrozenOnly
from tools.linear_algebra.gathering.regular.chain_matrix.main import Gathering_Matrix, Gathering_Vector



class _3dCSCG_TrForm_Numbering_Naive(FrozenOnly):
    def __init__(self, Tr):
        self._Tr_ = Tr
        self._mesh_ = Tr.mesh
        self._freeze_self_()

    def _3dCSCG_0Tr(self):
        """Do the numbering if it is a Tr 0-form

        :returns: A tuple of 4 outputs:

            1. (Gathering_Matrix) -- The global numbering in local mesh elements.
            2. dict -- The global numbering in local trace elements.
            3. (int) -- Number of dofs in this core.
            4. (None,...) -- Extra numbering information.
        """
        if self._Tr_.numbering._parameters_ == dict():
            return self._0Tr_no_parameters()
        else:
            raise NotImplementedError()


    def _0Tr_no_parameters(self):
        """
        Do the numbering if it is a Tr 0-form

        :returns: A tuple of 4 outputs:

            1. (Gathering_Matrix) -- The global numbering in local mesh elements.
            2. dict -- The global numbering in local trace elements.
            3. (int) -- Number of dofs in this core.
            4. (None,...) -- Extra numbering information.
        """
        GM = dict()
        GM_TEW = dict()
        local_num_dofs = 0
        extraInfo = None

        return GM, GM_TEW, local_num_dofs, extraInfo





    def _3dCSCG_1Tr(self):
        """Do the numbering if it is a Tr 1-form

        :returns: A tuple of 4 outputs:

            1. (Gathering_Matrix) -- The global numbering in local mesh elements.
            2. dict -- The global numbering in local trace elements.
            3. (int) -- Number of dofs in this core.
            4. (None,...) -- Extra numbering information.
        """
        if self._Tr_.numbering._parameters_ == dict():
            return self._1Tr_no_parameters()
        else:
            raise NotImplementedError()


    def _1Tr_no_parameters(self):
        """
        Do the numbering if it is a Tr 1-form

        :returns: A tuple of 4 outputs:

            1. (Gathering_Matrix) -- The global numbering in local mesh elements.
            2. dict -- The global numbering in local trace elements.
            3. (int) -- Number of dofs in this core.
            4. (None,...) -- Extra numbering information.
        """
        GM = dict()
        GM_TEW = dict()
        local_num_dofs = 0
        extraInfo = None

        return GM, GM_TEW, local_num_dofs, extraInfo









    def _3dCSCG_2Tr(self):
        """Do the numbering if it is a Tr 2-form

        :returns: A tuple of 4 outputs:

            1. (Gathering_Matrix) -- The global numbering in local mesh elements.
            2. dict -- The global numbering in local trace elements.
            3. (int) -- Number of dofs in this core.
            4. (None,...) -- Extra numbering information.
        """
        if self._Tr_.numbering._parameters_ == dict():
            return self._2Tr_no_parameters()
        else:
            raise NotImplementedError()


    def _2Tr_no_parameters(self):
        """
        Do the numbering if it is a Tr 2-form.

        :returns: A tuple of 4 outputs:

            1. (Gathering_Matrix) -- The global numbering in terms of local mesh elements.
            2. dict -- The global numbering in terms of local trace elements.
            3. (int) -- Number of dofs in this core.
            4. (None,...) -- Extra numbering information.
        """
        GM = dict()
        GM_TEW = dict()
        local_num_dofs = 0
        extraInfo = None

        num_basis_onside = self._Tr_.num.basis_onside
        NBO = [num_basis_onside['N'], num_basis_onside['W'], num_basis_onside['B']]

        type_amount_dict = self._mesh_.trace.elements.___PRIVATE_find_type_and_amount_numbered_before___()

        for i in self._mesh_.trace.elements:
            t_e_i = self._mesh_.trace.elements[i]
            am_NS, am_WE, am_BF = type_amount_dict[i]
            start_num = am_NS * NBO[0] + am_WE * NBO[1] +  am_BF * NBO[2]
            GM_TEW[i] = range(start_num, start_num + num_basis_onside[t_e_i.CHARACTERISTIC_side])
            local_num_dofs += num_basis_onside[t_e_i.CHARACTERISTIC_side]

        MAP = self._mesh_.trace.elements.map
        for i in MAP:
            vector = tuple()
            for t_j in MAP[i]:
                vector += (GM_TEW[t_j],)
            GM[i] = Gathering_Vector(i, vector)
        GM = Gathering_Matrix(GM, mesh_type='_3dCSCG')

        for i in GM_TEW:
            GM_TEW[i] = Gathering_Vector(i, GM_TEW[i])

        return GM, GM_TEW, local_num_dofs, extraInfo