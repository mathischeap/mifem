
from screws.frozen import FrozenOnly
from _3dCSCG.forms.standard._1_form.project.to import ___3dCSCG_1Form_Project_To___
from _3dCSCG.forms.standard._1_form.project.van import ___3dCSCG_1Form_Project_Van___



class _1Form_Projection(FrozenOnly):
    """A wrapper of all projection methods."""
    def __init__(self,_1sf):
        self._sf_ = _1sf
        self._to_ = ___3dCSCG_1Form_Project_To___(self._sf_)
        self._van_ = ___3dCSCG_1Form_Project_Van___(self._sf_)
        self._freeze_self_()

    @property
    def to(self):
        return self._to_

    @property
    def van(self):
        return self._van_


