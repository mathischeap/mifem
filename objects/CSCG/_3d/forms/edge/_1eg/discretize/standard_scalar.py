
from screws.freeze.base import FrozenOnly


class _3dCSCG_Edge1Form_Discretize_StandardScalar(FrozenOnly):
    """"""
    def __init__(self, ef):
        """"""
        self._ef_ = ef
        self._freeze_self_()


    def __call__(self, update_cochain=True, target='func'):
        """Discretize the standard _3dCSCG_ScalarField to a 1-edge-form

        'locally full local EEW cochain' means the cochain is a dict whose keys are edge-element
        numbers and values are edge-element-wise local cochains.

        """
        raise NotImplementedError()