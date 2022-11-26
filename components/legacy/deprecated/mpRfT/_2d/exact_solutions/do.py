# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/17 4:37 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly
from root.config.main import COMM, MPI, np


class mpRfT2_ES_Do(FrozenOnly):
    """"""

    def __init__(self, es):
        """"""
        self._es_ = es
        self._freeze_self_()


    def compute_Ln_norm_of(self, what, time=0, n=2, dp=2):
        """We compute the :math:`L^n`-norm of the attribute name ``what`` at `time`.

        :param str what:
        :param float time:
        :param int n: We compute :math:`L^n`-norm.
        :param dp:
        :return:
        """
        mesh = self._es_._mesh_
        what = getattr(self._es_.status, what)
        what.current_time = time

        assert mesh is not None, " <MS> : to compute L2_norm, I need a mesh."

        coo_map = mesh.coo_map.Gauss(dp)

        _, v = what.reconstruction(coo_map)
        detJ = mesh.rcMC.Jacobian(coo_map)

        LOCAL = 0
        for rc_rp in v:
            quad_weights = coo_map[rc_rp][1][1]
            local = np.sum([vij**n for vij in v[rc_rp]], axis=0)
            LOCAL += np.einsum('ij, ij, ij -> ', local, detJ[rc_rp], quad_weights, optimize='optimal')

        GLOBAL = COMM.allreduce(LOCAL, op=MPI.SUM) ** (1 / n)

        return GLOBAL


    def generate_random_valid_time_instances(self, amount=None):
        """

        Parameters
        ----------
        amount : None

        Returns
        -------

        """
        return self._es_.status.___Pr_generate_random_valid_time_instances___(
            amount=amount)


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
