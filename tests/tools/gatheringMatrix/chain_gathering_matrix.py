# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/29 2:26 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from components.freeze.base import FrozenOnly
from components.miscellaneous.miprint import miprint
from components.miscellaneous.mirand import sample
from tests.objects.CSCG._3d.randObj.form_caller import random_FormCaller_of_total_load_around as rf3
from objects.miUsGrid.triangular.master import Call as miUsFc2

from tools.linearAlgebra.gathering.regular.chain_matrix.main import Chain_Gathering_Matrix

class Test_ChainGM_sequent_chain_method(FrozenOnly):
    """"""

    def __init__(self):
        """"""
        miprint("SequentChainMethod [Test_ChainGM_sequent_chain_method] ...... ", flush=True)
        self.fc3 = rf3(500)
        self.fc2 = miUsFc2('st8', 2)
        self._freeze_self_()

    def __call__(self):
        """"""
        mesh = self.fc2.mesh
        f0 = self.fc2('0-f-i')
        f1 = self.fc2('1-f-i')
        f2 = self.fc2('2-f-o')
        gm0 = f0.numbering.gathering
        gm1 = f1.numbering.gathering
        gm2 = f2.numbering.gathering
        cgm = Chain_Gathering_Matrix([gm0, gm1, gm2], chain_method='sequent')
        check_dict_f0 = dict()
        check_dict_f1 = dict()
        check_dict_f2 = dict()
        for i in mesh.elements:
            fv = gm0[i].full_vector
            for j, sep_numbering in enumerate(fv):
                total_numbering = cgm[i][j]
                if sep_numbering not in check_dict_f0:
                    check_dict_f0[sep_numbering] = total_numbering
                else:
                    assert check_dict_f0[sep_numbering] == total_numbering

            fv = gm1[i].full_vector
            for j, sep_numbering in enumerate(fv):
                total_numbering = cgm[i][j + f0.num.basis]
                if sep_numbering not in check_dict_f1:
                    check_dict_f1[sep_numbering] = total_numbering
                else:
                    assert check_dict_f1[sep_numbering] == total_numbering

            fv = gm2[i].full_vector
            for j, sep_numbering in enumerate(fv):
                total_numbering = cgm[i][j + f0.num.basis + f1.num.basis]
                if sep_numbering not in check_dict_f2:
                    check_dict_f2[sep_numbering] = total_numbering
                else:
                    raise Exception(f"We should never reach this place.")

        #---- test find ----------------------------------------------------------------------------
        GLOBAL_num_dofs = cgm.GLOBAL_num_dofs
        sample_amount = int(GLOBAL_num_dofs/5)
        if sample_amount > 10: sample_amount = 10
        samples = sample(range(GLOBAL_num_dofs), sample_amount)

        for dof in samples:
            elements_indices = cgm.do.find.elements_and_local_indices_of_dof(dof)

            if elements_indices is not None:
                elements, indices = elements_indices

                for i in mesh.elements:
                    fv = cgm[i]

                    if i in elements:
                        index = elements.index(i)
                        local_numbering = indices[index]

                        assert fv[local_numbering] == dof
                    else:
                        assert dof not in fv
            else: # no local element found

                for i in mesh.elements:
                    fv = cgm[i]

                    assert dof not in fv

        elements_, indices_ = cgm.do.find.elements_and_local_indices_of_dofs(samples)
        for dof in samples:
            elements, indices = elements_[dof], indices_[dof]

            for i in mesh.elements:
                fv = cgm[i]

                if i in elements:
                    index = elements.index(i)
                    local_numbering = indices[index]
                    assert fv[local_numbering] == dof
                else:
                    assert dof not in fv

        #-------------------------------------------------------------------------------------------
        mesh = self.fc3.mesh
        f0 = self.fc3('0-f', is_hybrid=False)
        f1 = self.fc3('1-f', is_hybrid=False)
        f2 = self.fc3('2-f', is_hybrid=False)
        gm0 = f0.numbering.gathering
        gm1 = f1.numbering.gathering
        gm2 = f2.numbering.gathering
        cgm = Chain_Gathering_Matrix([gm0, gm1, gm2], chain_method='sequent')
        check_dict_f0 = dict()
        check_dict_f1 = dict()
        check_dict_f2 = dict()
        for i in mesh.elements:
            fv = gm0[i].full_vector
            for j, sep_numbering in enumerate(fv):
                total_numbering = cgm[i][j]
                if sep_numbering not in check_dict_f0:
                    check_dict_f0[sep_numbering] = total_numbering
                else:
                    assert check_dict_f0[sep_numbering] == total_numbering

            fv = gm1[i].full_vector
            for j, sep_numbering in enumerate(fv):
                total_numbering = cgm[i][j + f0.num.basis]
                if sep_numbering not in check_dict_f1:
                    check_dict_f1[sep_numbering] = total_numbering
                else:
                    assert check_dict_f1[sep_numbering] == total_numbering

            fv = gm2[i].full_vector
            for j, sep_numbering in enumerate(fv):
                total_numbering = cgm[i][j + f0.num.basis + f1.num.basis]
                if sep_numbering not in check_dict_f2:
                    check_dict_f2[sep_numbering] = total_numbering
                else:
                    assert check_dict_f2[sep_numbering] == total_numbering

        #---- test find ----------------------------------------------------------------------------
        GLOBAL_num_dofs = cgm.GLOBAL_num_dofs
        sample_amount = int(GLOBAL_num_dofs/5)
        if sample_amount > 10: sample_amount = 10
        samples = sample(range(GLOBAL_num_dofs), sample_amount)

        for dof in samples:
            elements_indices = cgm.do.find.elements_and_local_indices_of_dof(dof)

            if elements_indices is not None:
                elements, indices = elements_indices

                for i in mesh.elements:
                    fv = cgm[i]

                    if i in elements:
                        index = elements.index(i)
                        local_numbering = indices[index]

                        assert fv[local_numbering] == dof
                    else:
                        assert dof not in fv
            else: # no local element found

                for i in mesh.elements:
                    fv = cgm[i]

                    assert dof not in fv

        elements_, indices_ = cgm.do.find.elements_and_local_indices_of_dofs(samples)
        for dof in samples:
            elements, indices = elements_[dof], indices_[dof]

            for i in mesh.elements:
                fv = cgm[i]

                if i in elements:
                    index = elements.index(i)
                    local_numbering = indices[index]
                    assert fv[local_numbering] == dof
                else:
                    assert dof not in fv

        return 1


if __name__ == "__main__":
    # mpiexec -n 4 python __tests__/unittests/gathering_matrix/chain_gathering_matrix.py
    Test_ChainGM_sequent_chain_method()()
