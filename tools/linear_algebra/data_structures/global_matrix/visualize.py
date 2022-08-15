# -*- coding: utf-8 -*-




from screws.freeze.main import FrozenOnly
from root.config.main import rAnk, mAster_rank
import matplotlib.pyplot as plt




class ___GM_VISUALIZE___(FrozenOnly):
    def __init__(self, gm):
        self._gm_ = gm
        self._freeze_self_()

    def spy(self, markerfacecolor='k', markeredgecolor='g', markersize=6):
        """
        The spy plot of self.

        :param markerfacecolor:
        :param markeredgecolor:
        :param markersize:
        :return:
        """
        M = self._gm_.do.gather_M_to_core(core=mAster_rank)
        if rAnk == mAster_rank:
            fig = plt.figure()
            plt.spy(M,
                    markerfacecolor=markerfacecolor,
                    markeredgecolor=markeredgecolor,
                    markersize=markersize)
            plt.tick_params(axis='both', which='major', direction='out')
            plt.tick_params(which='both', top=True, right=True, labelbottom=True, labelright=True)
            plt.show()
            return fig
        else:
            return None

