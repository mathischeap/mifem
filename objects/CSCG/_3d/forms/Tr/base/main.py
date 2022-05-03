

from objects.CSCG._3d.forms.base import _3dCSCG_FORM_BASE
from objects.CSCG._3d.forms.Tr.base.IS import _3dCSCG_Tr_form_IS
from objects.CSCG._3d.forms.Tr.base.num import _3dCSCG_Tr_form_NUM
from objects.CSCG._3d.forms.Tr.base.numbering.main import _3dCSCG_Tr_Numbering
from objects.CSCG._3d.forms.Tr.base.cochain import _3dCSCG_Tr_Cochain
from objects.CSCG._3d.forms.Tr.base.do import _3dCSCG_Tr_form_DO
from objects.CSCG._3d.forms.Tr.base.matrices import _3dCSCG_Tr_Matrices





class _3dCSCG_Standard_Tr(_3dCSCG_FORM_BASE, ndim=3):
    """This is the parent of all 3d standard trace forms.

    :param mesh:
    :param space:
    :param str orientation: 'inner' or 'outer'.
    :param numbering_parameters: The parameters for the numbering. Including scheme name and other parameters.
        When it is a string, we use it as scheme name, and it has not other parameters.
    :type numbering_parameters: dict, str
    :param str name:
    """
    def __init__(self, mesh, space, orientation, numbering_parameters, name):
        """"""
        super().__init__(mesh, space)
        self._NUM_basis_, self._NUM_basis_components_, self._NUM_basis_onside_ = \
            getattr(self.space.num_basis, self.__class__.__name__)
        assert orientation in ('inner', 'outer'), " orientation needs to be 'inner' or 'outer'."
        self._orientation_ = orientation
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_Tr_form')
        self.standard_properties.name = name

        self._IS_ = None
        self._num_ = None
        self._numbering_ = _3dCSCG_Tr_Numbering(self, numbering_parameters)
        self._cochain_ = None
        self._do_ = None
        self._matrices_ = None


    def ___PRIVATE_reset_cache___(self):
        pass

    @property
    def orientation(self):
        """(str) The orientation."""
        return self._orientation_

    @property
    def IS(self):
        if self._IS_ is None:
            self._IS_ = _3dCSCG_Tr_form_IS(self)
        return self._IS_

    @property
    def num(self):
        if self._num_ is None:
            self._num_ = _3dCSCG_Tr_form_NUM(self)
        return self._num_

    @property
    def numbering(self):
        return self._numbering_

    @property
    def cochain(self):
        if self._cochain_ is None:
            self._cochain_ = _3dCSCG_Tr_Cochain(self)
        return self._cochain_

    @property
    def do(self):
        if self._do_ is None:
            self._do_ = _3dCSCG_Tr_form_DO(self)
        return self._do_

    @property
    def matrices(self):
        if self._matrices_ is None:
            self._matrices_ = _3dCSCG_Tr_Matrices(self)
        return self._matrices_