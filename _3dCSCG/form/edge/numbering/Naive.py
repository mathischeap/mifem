# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from SCREWS.frozen import FrozenOnly
from TOOLS.linear_algebra.gathering import Gathering_Matrix, Gathering_Vector



class _3dCSCG_Edge_Numbering_Naive(FrozenOnly):
    def __init__(self, ef):
        self._ef_ = ef
        self._mesh_ = ef.mesh
        self._freeze_self_()



    def _0Edge(self):
        """Do the numbering if it is an edge 0-form:
        :class:`_3dCSCG.form.edge._0_edge._0Edge`.

        :returns: A tuple of 4 outputs:

            1. (Gathering_Matrix) -- The global numbering in local elements.
            2. (Gathering_Matrix) -- The global numbering in local trace elements.
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
            2. (Gathering_Matrix) -- The global numbering in local edge elements. Edge-element-wise.
            3. (int) -- Number of dofs in this core.
            4. (None,...) -- Extra numbering information.
        """
        GM = dict()
        GM_EEW = dict()
        local_num_dofs = 0
        extraInfo = None

        # num_basis_onside = self._tf_.NUM_basis_onside
        # NBO = [num_basis_onside['N'], num_basis_onside['W'],  num_basis_onside['B']]
        #
        # type_amount_dict = self._mesh_.trace.elements.___DO_find_type_and_amount_numbered_before___()
        #
        # for i in self._mesh_.trace.elements:
        #     t_e_i = self._mesh_.trace.elements[i]
        #     am_NS, am_WE, am_BF = type_amount_dict[i]
        #     start_num = am_NS * NBO[0] + am_WE * NBO[1] +  am_BF * NBO[2]
        #     GM_TEW[i] = range(start_num, start_num + num_basis_onside[t_e_i.CHARACTERISTIC_side])
        #
        # MAP = self._mesh_.trace.elements.map
        # for i in MAP:
        #     vector = tuple()
        #     for t_j in MAP[i]:
        #         vector += (GM_TEW[t_j],)
        #     GM[i] = Gathering_Vector(i, vector)
        # GM = Gathering_Matrix(GM, mesh_type='_3dCSCG')
        #
        # for i in GM_TEW:
        #     GM_TEW[i] = Gathering_Vector(i, GM_TEW[i])

        return GM, GM_EEW, local_num_dofs, extraInfo