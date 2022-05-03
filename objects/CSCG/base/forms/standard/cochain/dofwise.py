# -*- coding: utf-8 -*-
from screws.freeze.base import FrozenOnly
import numpy as np
from root.config.main import cOmm, rAnk, mAster_rank



class CSCG_SF_Cochain_DofWise(FrozenOnly):
    """"""
    def __init__(self, cochain):
        """"""
        self._cochain_ = cochain
        self._f_ = cochain._sf_
        self._freeze_self_()

    def __getitem__(self, i):
        """Return the cochain of the global dof #i in all cores."""
        GM = self._f_.numbering.gathering
        assert i % 1 == 0 and -self._f_.num.GLOBAL_dofs <= i < self._f_.num.GLOBAL_dofs, f'i={i} is wrong!'

        if i < 0:
            i += self._f_.num.GLOBAL_dofs

        ME_LC = GM.do.find.elements_and_local_indices_of_dof(i)

        cc = None

        if ME_LC is None:
            pass
        else:
            mesh_elements, local_indices = ME_LC
            for e, ind in zip(mesh_elements, local_indices):
                if cc is None:
                    cc = self._cochain_.local[e][ind]
                else:
                    np.testing.assert_almost_equal(cc, self._cochain_.local[e][ind])

        CC = cOmm.gather(cc, root=mAster_rank)

        if rAnk == mAster_rank:
            the_cc = None
            for cc in CC:
                if cc is not None:
                    if the_cc is None:
                        the_cc = cc
                    else:
                        np.testing.assert_almost_equal(cc, the_cc)
                else:
                    pass
            assert the_cc is not None, f"We must have found a cochain."
        else:
            the_cc = None

        return cOmm.bcast(the_cc, root=mAster_rank)