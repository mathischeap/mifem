# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/28/2022 12:01 AM
"""
import sys, os

if './' not in sys.path: sys.path.append('./')
from screws.freeze.main import FrozenOnly


from playsound import playsound
from root.config.main import rAnk, mAster_rank



class MiSoundEffectError(Exception):
    """Raise when we try to define new attribute for a frozen object."""


class MiSoundEffect(FrozenOnly):
    """"""

    def __init__(self, ):
        """"""
        self._freeze_self_()


    @classmethod
    def transition(cls, i=0):
        """Play the ith transition sound effect.

        Parameters
        ----------
        i

        Returns
        -------

        """

        # noinspection PyBroadException
        if rAnk != mAster_rank: return

        try:
            absolute_path = os.path.dirname(__file__)

            if i == 0:
                playsound(absolute_path + '/mixkit_fast_small_sweep_transition_166.wav')
            elif i == 1:
                playsound(absolute_path + '/transitions/mixkit-arcade-retro-game-over-213.wav')
            elif i == 2:
                playsound(absolute_path + '/transitions/mixkit-cartoon-toy-whistle-616.wav')
            elif i == 3:
                playsound(absolute_path + '/transitions/mixkit-classic-alarm-995.wav')
            elif i == 4:
                playsound(absolute_path + '/transitions/mixkit-crowd-laugh-424.wav')
            elif i == 5:
                playsound(absolute_path + '/transitions/mixkit-dog-barking-twice-1.wav')
            elif i == 6:
                playsound(absolute_path + '/transitions/mixkit-fast-rocket-whoosh-1714.wav')
            elif i == 7:
                playsound(absolute_path + '/transitions/mixkit-retro-game-notification-212.wav')
            elif i == 8:
                playsound(absolute_path + '/transitions/mixkit-sad-game-over-trombone-471.wav')
            elif i == 9:
                playsound(absolute_path + '/transitions/mixkit-small-crowd-laugh-and-applause-422.wav')
            else:
                raise MiSoundEffectError(f"{i}th transition sound effect is not found.")

        except:
            pass



if __name__ == '__main__':
    # mpiexec -n 4 python screws/sound_effects/main.py
    MiSoundEffect.transition(0)
