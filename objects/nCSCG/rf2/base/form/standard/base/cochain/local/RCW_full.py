# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/16 6:00 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.nCSCG.rf2.base.form.standard.base.cochain.local.base import nCSCG_RF2_LocalCochainBase





class nCSCG_RF2__RCW_Full__LocalCochain(nCSCG_RF2_LocalCochainBase):
    """root-cell-wise; all root-cells are available."""

    def __init__(self, f, LC):
        """"""
        super(nCSCG_RF2__RCW_Full__LocalCochain, self).__init__(f, LC)
        self._freeze_self_()


    def ___Pr_check_and_parse_LC___(self, LC):
        """"""
        mesh = self._mesh_
        for i in mesh: assert repr(mesh(i)) in LC
        return LC

    def __getitem__(self, indices):
        """

        Parameters
        ----------
        indices :
            the indices of a root cell.

        Returns
        -------

        """
        assert self.signature == self._mesh_.signature, f"signature dis-match."
        return self._LC_[repr(self._mesh_(indices))]




if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
