
from screws.freeze.main import FrozenOnly
from objects.CSCG._3d.forms.standard._2s.project.to import ___3dCSCG_2Form_Project_To___
from objects.CSCG._3d.forms.standard._2s.project.van import ___3dCSCG_2Form_Project_Van___


class _2Form_Projection(FrozenOnly):
    """A wrapper of all projection methods."""
    def __init__(self,_2sf):
        self._sf_ = _2sf
        self._to_ = ___3dCSCG_2Form_Project_To___(self._sf_)
        self._van_ = ___3dCSCG_2Form_Project_Van___(self._sf_)
        self._freeze_self_()

    @property
    def to(self):
        return self._to_

    @property
    def van(self):
        return self._van_