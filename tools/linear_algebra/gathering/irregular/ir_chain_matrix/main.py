# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/23 10:27 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from root.config.main import cOmm, MPI, np
from tools.linear_algebra.gathering.irregular.ir_chain_matrix.do.main import iR_CGM_DO



class iR_Chain_Gathering_Matrix(FrozenOnly):
    """"""

    def __init__(self, GMs):
        """
        GMs can be a gathering matrix or a chain gathering matrix. Or a list of them.

        :param GMs:
        """

        if GMs.__class__.__name__ in ('iR_Gathering_Matrix', 'iR_Chain_Gathering_Matrix,'):
            GMs = [GMs,]

        assert isinstance(GMs, (list, tuple)) and len(GMs) >= 1, \
            "we need a list or tuple of at least one irregular gathering matrices."

        mesh_types = list()
        NEW_GMs = list()
        for i, gm in enumerate(GMs):

            mesh_types.append(gm.mesh_type)
            assert mesh_types[i] == mesh_types[0], \
                f"mesh_types[{i}]: {mesh_types[i]} is different from mesh_types[0]: {mesh_types[0]}"

            if gm.__class__.__name__ == 'iR_Gathering_Matrix':
                NEW_GMs.append(gm)
            elif gm.__class__.__name__ == 'iR_Chain_Gathering_Matrix':
                NEW_GMs.extend(gm.GMs)
            else:
                raise Exception(f'GMs[{i}] is {gm.__class__.__name__}, wrong!')

        self._mesh_type_ = mesh_types[0]
        GMs = NEW_GMs

        self.___NUM___ = len(GMs)
        # ... check GM in GMs are actually representing same local elements.
        if self.___NUM___ == 1:
            pass
        else:
            gm0 = GMs[0]
            LEN = len(gm0)
            for _, gm in enumerate(GMs[1:]):
                assert len(gm) == LEN, f"GMs[{_+1}] has different length comparing to GMs[0]."
                for i in gm0:
                    assert i in gm, f"GMs[{_+1}] represent different elements comparing to GMs[0]."

        To_Be_Added = [0,]
        for gm in GMs[:-1]:
            GLOBAL_num_dofs = gm.GLOBAL_num_dofs
            # noinspection PyUnresolvedReferences
            To_Be_Added.append(To_Be_Added[-1]+GLOBAL_num_dofs)

        self._To_Be_Added_ = To_Be_Added
        self._GLOBAL_num_dofs_ = To_Be_Added[-1] + GMs[-1].GLOBAL_num_dofs
        self._GMs_ = GMs
        self._NUM_GMs_ = len(GMs)
        self._dofs_distribution_ = [_.GLOBAL_num_dofs for _ in self.GMs]

        self.___LENs_pool___ = list()
        self.___LENs___ = dict()
        self._local_ranges_ = None
        self._do_ = iR_CGM_DO(self)
        self._freeze_self_()

    @property
    def ___Pr_IS_regular___(self):
        return False

    @property
    def NUM_GMs(self):
        """Return an int which represents how many GM this CGM has."""
        return self._NUM_GMs_

    @property
    def GLOBAL_num_dofs(self):
        """How many dofs in total in all cores."""
        return self._GLOBAL_num_dofs_


    def __getitem__(self, item):
        """Return the chained gathering vector for the mesh element numbered by `i`.

        :param item: #i mesh element.
        :return: A 1d numpy.array. But we do not save it, it is made whenever we call this method.
        """
        if self.___NUM___ == 1:
            return self.GMs[0][item].full_vector
        else:
            _ = [gm[item].full_vector+self._To_Be_Added_[i] for i, gm in enumerate(self._GMs_)]
            return np.concatenate(_)

    def __iter__(self):
        """Go through all local mesh elements."""
        for rp in self._GMs_[0]:
            yield rp

    def __contains__(self, rp):
        """If a mesh element is contained by this CGM? Or if a mesh element is local?"""
        return rp in self._GMs_[0]

    def __len__(self):
        """local length, equal to the amount of local elements."""
        return len(self._GMs_[0])

    def __eq__(self, other):
        """If two CGMs are equal?"""
        if self is other:
            RETURN = True
        else:
            if other.__class__.__name__ != 'iR_Chain_Gathering_Matrix':
                RETURN = False
            elif other.___NUM___ != self.___NUM___:
                RETURN = False
            else:
                RETURN = True
                for i, gm in enumerate(self._GMs_):
                    tf = gm == other._GMs_[i]
                    if not tf:
                        RETURN = False
                        break

        RETURN = cOmm.allreduce(RETURN, op=MPI.LAND)
        return RETURN

    @property
    def mesh_type(self):
        return self._mesh_type_

    @property
    def GMs(self):
        return self._GMs_

    def ___Pr_LENs___(self, rc_rp):
        """Return the length of each local GM in the root-cell."""
        if rc_rp not in self.___LENs___:

            LENs = list()
            for gm in self.GMs:
                LENs.append(len(gm[rc_rp]))
            if LENs in self.___LENs_pool___:
                i = self.___LENs_pool___.index(LENs)
                self.___LENs___[rc_rp] = self.___LENs_pool___[i]
            else:
                self.___LENs_pool___.append(LENs)
                self.___LENs___[rc_rp] = self.___LENs_pool___[-1]

        return self.___LENs___[rc_rp]

    @property
    def local_ranges(self):
        """The local dofs are in this range."""
        if self._local_ranges_ is None:

            if len(self) == 0:
                self._local_ranges_ = [tuple() for _ in range(len(self.GMs))]

            else:
                self._local_ranges_ = list()
                tba = self._To_Be_Added_
                for i, gm in enumerate(self.GMs):
                    LR = gm.local_range
                    MIN, MAX = LR
                    MIN += tba[i]
                    MAX += tba[i]
                    self._local_ranges_.append((MIN, MAX))

        return self._local_ranges_

    @property
    def GLOBAL_dofs_distribution(self):
        """Like [100, 127, 80], we know GMs[0].GLOBAL_num_dofs=100, GMs[1].GLOBAL_num_dofs=127,
         GMs[2].GLOBAL_num_dofs=80. So we know, for example, how to distribution the result vector
         by solving the assembled EWC_SparseMatrix."""
        return self._dofs_distribution_

    @property
    def do(self):
        return self._do_






if __name__ == "__main__":
    # mpiexec -n 4 python tools/linear_algebra/gathering/irregular/ir_chain_matrix/main.py
    from __init__ import rfT2
    fc = rfT2.rf(100)

    f2 = fc('2-f-o')
    f1 = fc('1-f-o')
    nt = fc('nst')

    GM2 = f2.numbering.gathering
    GM1 = f1.numbering.gathering
    GMT = nt.numbering.gathering

    GM1 = iR_Chain_Gathering_Matrix([GM2, GM1, GMT])
    GM2 = iR_Chain_Gathering_Matrix(GM2)


