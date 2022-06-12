# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/31 9:21 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.mpRfT._2d.mesh.space.helpers.base import mpRfT2_Mesh_Space_Basis





class mpRfT2_Mesh_Space_NSgF_Basis(mpRfT2_Mesh_Space_Basis):
    def __init__(self, t, mesh, coo_map):
        """"""
        super(mpRfT2_Mesh_Space_NSgF_Basis, self).__init__(mesh, coo_map)
        self._t_ = t
        self._freeze_self_()

    @property
    def ___Pr_rcMC_key___(self):
        """Cannot be used for rcMC"""
        raise Exception("Cannot be used for rcMC")

    def ___Pr_rcMC_nodes___(self, rp):
        """Cannot be used for rcMC"""
        raise Exception("Cannot be used for rcMC")

    @property
    def ___Pr_sgMC_key___(self):
        """"""
        raise NotImplementedError()

    def ___Pr_sgMC_nodes___(self, rp):
        """"""
        raise NotImplementedError()

    def __getitem__(self, seg):
        """Get the basis functions for the root-cell mesh(indices)"""
        assert seg.__class__.__name__ == 'mpRfT2_Segment', f"I need a mpRfT2 segment."
        return self._getitem_(seg)

    def ___Pr_getitem_uniform___(self, seg):
        """"""
        rp = seg.__repr__()
        N = self._t_.N[seg]

        if rp[3] == 'c':
            ind = rp.split(':')[-1]
            N_ind = str(N) + ind
        elif rp[3] == 't':
            ind =  rp.split('-')[-1]
            N_ind = str(N) + ind
        else:
            raise Exception()

        if N_ind in self._cache_:
            pass
        else:
            space = self._mesh_.space[N]
            nodes = self._cm_[seg]

            node_basis = space.basises[0].node_basis(x=nodes)

            _basis_ = (node_basis,)

            self._cache_[N_ind] = nodes, _basis_

        return self._cache_[N_ind]

    def ___Pr_getitem_Gauss___(self, seg):
        return self.___Pr_getitem_Lobatto___(seg)

    def ___Pr_getitem_Lobatto___(self, seg):
        N = self._t_.N[seg]

        if N in self._cache_:
            pass
        else:
            space = self._mesh_.space[N]
            nodes = self._cm_[seg][0] # get nodes, [1] is the weights.

            node_basis = space.basises[0].node_basis(x=nodes)

            _basis_ = (node_basis,)

            self._cache_[N] = nodes, _basis_

        return self._cache_[N]







if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/mesh/space/helpers/nsg.py
    pass
