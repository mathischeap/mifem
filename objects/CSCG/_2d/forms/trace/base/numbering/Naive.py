# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, the Netherlands

"""
from components.freeze.main import FrozenOnly
from tools.linearAlgebra.gathering.regular.chain_matrix.main import Gathering_Matrix, Gathering_Vector



class _2dCSCG_Trace_Numbering_Naive(FrozenOnly):
    def __init__(self, tf):
        self._tf_ = tf
        self._mesh_ = tf.mesh
        self._freeze_self_()

    def _2dCSCG_1Trace_Outer(self):
        """
        Do the numbering if it is a outer trace 1-form.

        :returns: A tuple of 3 outputs:

            1. (Gathering_Matrix) -- The global numbering in local elements.
            2. (Gathering_Matrix) -- The global numbering in local trace elements.
            3. (int) -- Number of dofs in this core.
            4. (None,...) -- Extra numbering information.
        """
        if self._tf_.numbering._parameters_ == dict():
            return self._1Trace_Outer_no_parameters()
        else:
            raise NotImplementedError()



    def _1Trace_Outer_no_parameters(self):
        """
        Do the numbering if it is a outer trace 1-form.

        :returns: A tuple of 3 outputs:

            1. (Gathering_Matrix) -- The global numbering in local elements.
            2. (Gathering_Matrix) -- The global numbering in local trace elements.
            3. (int) -- Number of dofs in this core.
            4. (None,...) -- Extra numbering information.
        """
        GM = dict()
        GM_TEW = dict()
        local_num_dofs = 0
        extraInfo = None

        num_basis_onside = self._tf_.num.basis_onside
        NBO = [num_basis_onside['U'], num_basis_onside['L']]

        type_amount_dict = self._mesh_.trace.elements.___PRIVATE_find_type_and_amount_numbered_before___()

        for i in self._mesh_.trace.elements:
            t_e_i = self._mesh_.trace.elements[i]
            am_UD, am_LR = type_amount_dict[i]
            start_num = am_UD * NBO[0] + am_LR * NBO[1]
            GM_TEW[i] = range(start_num, start_num + num_basis_onside[t_e_i.CHARACTERISTIC_edge])

        MAP = self._mesh_.trace.elements.map
        for i in MAP:
            vector = tuple()
            for t_j in MAP[i]:
                vector += (GM_TEW[t_j],)
            GM[i] = Gathering_Vector(i, vector)
        GM = Gathering_Matrix(GM, mesh_type='_2dCSCG')

        for i in GM_TEW:
            GM_TEW[i] = Gathering_Vector(i, GM_TEW[i])

        return GM, GM_TEW, local_num_dofs, extraInfo