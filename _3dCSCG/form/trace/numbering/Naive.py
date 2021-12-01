# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from SCREWS.frozen import FrozenOnly
from TOOLS.linear_algebra.gathering import Gathering_Matrix, Gathering_Vector



class _3dCSCG_Trace_Form_Numbering_Naive(FrozenOnly):
    def __init__(self, tf):
        self._tf_ = tf
        self._mesh_ = tf.mesh
        self._freeze_self_()

    def _0Trace(self):
        """Do the numbering if it is a trace 0-form:
        :class:`_3dCSCG.form.standard._0_trace._0Trace`.

        :returns: A tuple of 4 outputs:

            1. (Gathering_Matrix) -- The global numbering in local elements.
            2. dict -- The global numbering in local trace elements.
            3. (int) -- Number of dofs in this core.
            4. (None,...) -- Extra numbering information.
        """
        if self._tf_.numbering._parameters_ == dict():
            return self._0Trace_no_parameters()
        else:
            raise NotImplementedError()


    def _0Trace_no_parameters(self):
        """
        Do the numbering if it is a trace 0-form:
        :class:`_3dCSCG.form.standard._0_trace._0Trace`.

        :returns: A tuple of 4 outputs:

            1. (Gathering_Matrix) -- The global numbering in local elements.
            2. dict -- The global numbering in local trace elements.
            3. (int) -- Number of dofs in this core.
            4. (None,...) -- Extra numbering information.
        """
        GM = dict()
        GM_TEW = dict()
        local_num_dofs = 0
        extraInfo = None

        num_basis_onside = self._tf_.NUM_basis_onside
        NBO = [num_basis_onside['N'], num_basis_onside['W'],  num_basis_onside['B']]

        type_amount_dict = self._mesh_.trace.elements.___DO_find_type_and_amount_numbered_before___()

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





    def _1Trace(self):
        """Do the numbering if it is a trace 1-form:
        :class:`_3dCSCG.form.standard._1_trace._1Trace`.

        :returns: A tuple of 4 outputs:

            1. (Gathering_Matrix) -- The global numbering in local elements.
            2. dict -- The global numbering in local trace elements.
            3. (int) -- Number of dofs in this core.
            4. (None,...) -- Extra numbering information.
        """
        if self._tf_.numbering._parameters_ == dict():
            return self._0Trace_no_parameters()
        else:
            raise NotImplementedError()


    def _1Trace_no_parameters(self):
        """
        Do the numbering if it is a trace 1-form:
        :class:`_3dCSCG.form.standard._1_trace._1Trace`.

        :returns: A tuple of 4 outputs:

            1. (Gathering_Matrix) -- The global numbering in local elements.
            2. dict -- The global numbering in local trace elements.
            3. (int) -- Number of dofs in this core.
            4. (None,...) -- Extra numbering information.
        """
        GM = dict()
        GM_TEW = dict()
        local_num_dofs = 0
        extraInfo = None

        num_basis_onside = self._tf_.NUM_basis_onside
        NBO = [num_basis_onside['N'], num_basis_onside['W'],  num_basis_onside['B']]

        type_amount_dict = self._mesh_.trace.elements.___DO_find_type_and_amount_numbered_before___()

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









    def _2Trace(self):
        """Do the numbering if it is a trace 2-form:
        :class:`_3dCSCG.form.standard._2_trace._2Trace`.

        :returns: A tuple of 4 outputs:

            1. (Gathering_Matrix) -- The global numbering in local elements.
            2. dict -- The global numbering in local trace elements.
            3. (int) -- Number of dofs in this core.
            4. (None,...) -- Extra numbering information.
        """
        if self._tf_.numbering._parameters_ == dict():
            return self._2Trace_no_parameters()
        else:
            raise NotImplementedError()


    def _2Trace_no_parameters(self):
        """
        Do the numbering if it is a trace 2-form:
        :class:`_3dCSCG.form.standard._2_trace._2Trace`.

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

        num_basis_onside = self._tf_.NUM_basis_onside
        NBO = [num_basis_onside['N'], num_basis_onside['W'],  num_basis_onside['B']]

        type_amount_dict = self._mesh_.trace.elements.___DO_find_type_and_amount_numbered_before___()

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