from screws.freeze.base import FrozenOnly

from objects.CSCG._2d.forms.standard.base.dofs.dof.visualize import _2dCSCG_SF_DOF_VIS

class _2dCSCG_SF_DOF(FrozenOnly):
    """"""
    def __init__(self, dofs, i):
        """ I am the #i global dof of a standard 2d CSCG form.

        :param dofs:
        :param i:
        """

        self._dofs_ = dofs
        self._i_ = i
        self._visualize_ = None
        self._freeze_self_()

    @property
    def i(self):
        return self._i_

    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = _2dCSCG_SF_DOF_VIS(self)
        return self._visualize_

    @property
    def cochain(self):
        return self._dofs_._sf_.cochain.dofwise[self.i]

