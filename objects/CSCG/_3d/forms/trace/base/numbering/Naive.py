# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import numpy as np
from root.config.main import RANK, SIZE, COMM, MPI, MASTER_RANK
from components.freeze.main import FrozenOnly
from tools.elementwiseCache.gathering.regular.chain_matrix.main import Gathering_Matrix, Gathering_Vector


class _3dCSCG_Trace_Form_Numbering_Naive(FrozenOnly):
    """"""

    def __init__(self, tf):
        """"""
        self._tf_ = tf
        self._mesh_ = tf.mesh
        self._freeze_self_()

    def _3dCSCG_0Trace(self):
        """Do the numbering if it is a trace 0-form:
        :class:`_3dCSCG.form.standard._0_trace._0Trace`.

        :returns: A tuple of 4 outputs:

            1. (Gathering_Matrix) -- The global numbering in local elements.
            2. dict -- The global numbering in local trace elements.
            3. (int) -- Number of dofs in this core.
            4. (None,...) -- Extra numbering information.
        """
        if self._tf_.numbering._parameters_ == dict():
            return self._0Trace_no_parameters(self._tf_.whether.hybrid)
        else:
            raise NotImplementedError()

    def _0Trace_no_parameters(self, hybrid):
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
        GM_TEW: dict = dict()
        local_num_dofs = 0
        extraInfo = None

        if hybrid:
            num_basis_onside = self._tf_.num.basis_onside
            NBO = [num_basis_onside['N'], num_basis_onside['W'], num_basis_onside['B']]

            type_amount_dict = self._mesh_.trace.elements.\
                ___PRIVATE_find_type_and_amount_numbered_before___()

            for i in self._mesh_.trace.elements:
                t_e_i = self._mesh_.trace.elements[i]
                am_NS, am_WE, am_BF = type_amount_dict[i]
                start_num = am_NS * NBO[0] + am_WE * NBO[1] + am_BF * NBO[2]
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

        else:

            GM, GM_TEW, local_num_dofs, extraInfo = self.___Pr_nonhybrid_numbering___(0)

        return GM, GM_TEW, local_num_dofs, extraInfo

    def _3dCSCG_1Trace(self):
        """Do the numbering if it is a trace 1-form:
        :class:`_3dCSCG.form.standard._1_trace._1Trace`.

        :returns: A tuple of 4 outputs:

            1. (Gathering_Matrix) -- The global numbering in local elements.
            2. dict -- The global numbering in local trace elements.
            3. (int) -- Number of dofs in this core.
            4. (None,...) -- Extra numbering information.
        """
        if self._tf_.numbering._parameters_ == dict():
            return self._1Trace_no_parameters(self._tf_.whether.hybrid)
        else:
            raise NotImplementedError()


    def _1Trace_no_parameters(self, hybrid):
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
        GM_TEW: dict = dict()
        local_num_dofs = 0
        extraInfo = None

        if hybrid:
            num_basis_onside = self._tf_.num.basis_onside
            NBO = [num_basis_onside['N'], num_basis_onside['W'],  num_basis_onside['B']]

            type_amount_dict = self._mesh_.trace.elements.___PRIVATE_find_type_and_amount_numbered_before___()

            for i in self._mesh_.trace.elements:
                t_e_i = self._mesh_.trace.elements[i]
                am_NS, am_WE, am_BF = type_amount_dict[i]
                start_num = am_NS * NBO[0] + am_WE * NBO[1] + am_BF * NBO[2]
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

        else:

            GM, GM_TEW, local_num_dofs, extraInfo = self.___Pr_nonhybrid_numbering___(1)

        return GM, GM_TEW, local_num_dofs, extraInfo

    def _3dCSCG_2Trace(self):
        """Do the numbering if it is a trace 2-form:
        :class:`_3dCSCG.form.standard._2_trace._2Trace`.

        :returns: A tuple of 4 outputs:

            1. (Gathering_Matrix) -- The global numbering in local elements.
            2. dict -- The global numbering in local trace elements.
            3. (int) -- Number of dofs in this core.
            4. (None,...) -- Extra numbering information.
        """
        if self._tf_.numbering._parameters_ == dict():
            return self._2Trace_no_parameters(self._tf_.whether.hybrid)
        else:
            raise NotImplementedError()


    def _2Trace_no_parameters(self, hybrid):
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

        if hybrid:
            num_basis_onside = self._tf_.num.basis_onside
            NBO = [num_basis_onside['N'], num_basis_onside['W'], num_basis_onside['B']]

            type_amount_dict = self._mesh_.trace.elements.___PRIVATE_find_type_and_amount_numbered_before___()

            for i in self._mesh_.trace.elements:
                t_e_i = self._mesh_.trace.elements[i]
                am_NS, am_WE, am_BF = type_amount_dict[i]
                start_num = am_NS * NBO[0] + am_WE * NBO[1] + am_BF * NBO[2]
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

        else:

            GM, GM_TEW, local_num_dofs, extraInfo = self.___Pr_nonhybrid_numbering___(2)

        return GM, GM_TEW, local_num_dofs, extraInfo

    def ___Pr_nonhybrid_numbering___(self, k):
        """A general routine to get non-hybrid trace-form numbering from non-hybrid sf numbering."""

        GM = dict()
        GM_TEW = dict()
        extraInfo = None

        if k == 0:
            from objects.CSCG._3d.forms.standard._0s.main import _3dCSCG_0Form as form
            internal_local_dofs = self._tf_.space.internal_local_dofs._3dCSCG_0Form

        elif k == 1:
            from objects.CSCG._3d.forms.standard._1s.main import _3dCSCG_1Form as form
            internal_local_dofs = self._tf_.space.internal_local_dofs._3dCSCG_1Form

        elif k == 2:
            from objects.CSCG._3d.forms.standard._2s.main import _3dCSCG_2Form as form
            internal_local_dofs = self._tf_.space.internal_local_dofs._3dCSCG_2Form

        else:
            raise Exception()

        asf = form(self._tf_.mesh, self._tf_.space, numbering_parameters='Naive', hybrid=False)

        sfGM = asf.numbering.gathering
        global_num_dofs = sfGM.global_num_dofs
        num_basis = asf.num.basis

        mesh = self._mesh_
        local_GM = dict()
        for i in mesh.elements:
            local_GM[i] = sfGM[i].full_vector
        del sfGM

        amount_internal_dofs = len(internal_local_dofs)
        surface_local_dofs = list()
        for i in range(num_basis):
            if i in internal_local_dofs:
                pass
            else:
                surface_local_dofs.append(i)

        # we will let rank 0 to do all the numbering .........
        if RANK == 0:
            ReNumberingArray = np.zeros(global_num_dofs, dtype=int)
            CN = 0
            number2 = -1
            for i in local_GM:
                GV = local_GM[i]
                if amount_internal_dofs > 0:
                    ReNumberingArray[GV[internal_local_dofs]] = -1
                else:
                    pass

                surface_global_dofs = GV[surface_local_dofs]

                surface_global_dofs = surface_global_dofs[surface_global_dofs > number2]

                LEN = len(surface_global_dofs)
                if LEN == 0:
                    pass
                else:
                    ReNumberingArray[surface_global_dofs] = np.arange(CN, CN+LEN)
                    number2 = max(surface_global_dofs)
                    CN += LEN

            for rank in range(1, SIZE):

                local_GM_rank = COMM.recv(source=rank, tag=rank)

                for i in local_GM_rank:
                    GV = local_GM_rank[i]
                    if amount_internal_dofs > 0:
                        ReNumberingArray[GV[internal_local_dofs]] = -1
                    else:
                        pass

                    surface_global_dofs = GV[surface_local_dofs]

                    surface_global_dofs = surface_global_dofs[surface_global_dofs > number2]

                    LEN = len(surface_global_dofs)
                    if LEN == 0:
                        pass
                    else:
                        ReNumberingArray[surface_global_dofs] = np.arange(CN, CN+LEN)
                        number2 = max(surface_global_dofs)
                        CN += LEN

                COMM.Send([ReNumberingArray, MPI.INT], dest=rank, tag=rank+1)

            assert max(ReNumberingArray) + 1 == \
                   global_num_dofs - mesh.elements.global_num * amount_internal_dofs

        else:
            del internal_local_dofs

            COMM.send(local_GM, dest=0, tag=RANK)

            ReNumberingArray = np.zeros(global_num_dofs, dtype=int)
            COMM.Recv([ReNumberingArray, MPI.INT], source=0, tag=RANK+1)

        local_gathering = self._tf_.numbering.local_gathering

        from components.assemblers import VectorAssembler
        assembler = VectorAssembler(local_gathering)

        Tmap = mesh.trace.elements.map
        GM_temp = dict()

        all_local_dofs = set()

        for i in local_GM:

            tes = Tmap[i]

            for te, side in zip(tes, 'NSWEBF'):
                sf_side_local_dofs = asf.numbering.do.find.local_dofs_on_element_side(side)
                _ = ReNumberingArray[local_GM[i]][sf_side_local_dofs]
                GM_temp[side] = _
                GM_TEW[te] = Gathering_Vector(te, _)

                all_local_dofs.update(_)

            GM[i] = Gathering_Vector(i, assembler(GM_temp, 'replace', dtype=int))

        assert -1 not in all_local_dofs, f"internal dof found!"
        local_num_dofs = len(all_local_dofs)
        all_local_dofs = COMM.gather(all_local_dofs, root=MASTER_RANK)
        if RANK == MASTER_RANK:
            _ = set()
            for i in range(SIZE):
                _.update(all_local_dofs[i])
                all_local_dofs[i] = None
            assert min(_) == 0
        else:
            _ = None

        GM = Gathering_Matrix(GM, mesh_type='_3dCSCG')
        assert GM.global_num_dofs == global_num_dofs - mesh.elements.global_num * amount_internal_dofs
        if RANK == MASTER_RANK:
            MAX = max(_)
            assert MAX == GM.global_num_dofs - 1
            assert len(_) == MAX + 1

        return GM, GM_TEW, local_num_dofs, extraInfo
