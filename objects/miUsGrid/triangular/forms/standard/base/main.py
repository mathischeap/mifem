# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/20 2:36 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from objects.miUsGrid.triangular.forms.base.main import miUsTriangular_FormBase
from objects.miUsGrid.triangular.forms.standard.base.cochain.main import miUs_Triangular_SF_Cochain
from objects.miUsGrid.triangular.forms.standard.base.num import miUs_Triangular_SF_Num
from objects.miUsGrid.triangular.forms.standard.base.coboundary import miUs_Triangular_SF_Coboundary
from objects.miUsGrid.triangular.forms.standard.base.numbering.main import miUs_Triangular_SF_Numbering
from objects.miUsGrid.triangular.forms.standard.base.IS import miUs_Triangular_SF_IS
from objects.miUsGrid.triangular.forms.standard.base.BC.main import miUsTriangle_SF_BC


class miUsTriangular_SF_Base(miUsTriangular_FormBase):
    """"""

    def __init__(self, mesh, space, orientation, k, name):
        """"""
        super(miUsTriangular_SF_Base, self).__init__(mesh, space, name)
        assert orientation in ('outer', 'inner'), f"orientation={orientation} invalid."
        self.standard_properties.___PRIVATE_add_tag___('miUsGrid_triangular_standard_form')
        self._orientation_ = orientation
        assert k in (0, 1, 2), f"{k}-form is invalid."
        self._k_ = k

        self._cochain_ = miUs_Triangular_SF_Cochain(self)
        self._num_ = miUs_Triangular_SF_Num(self)
        self._coboundary_ = miUs_Triangular_SF_Coboundary(self)
        self._numbering_ = miUs_Triangular_SF_Numbering(self)
        self._IS_ = miUs_Triangular_SF_IS(self)
        self._BC_ = None

    def __repr__(self):
        """"""
        return f"miUsTriangular_S{self.k}F_{self.name}@{id(self)}"

    @property
    def shadow(self):
        # noinspection PyArgumentList
        shadow = self.__class__(self.mesh, self.space, name = 'shadow-of-' + self.name)
        shadow.numbering._gathering_ = self.numbering.gathering
        return  shadow

    @property
    def ___Pr_EWC_cache_key___(self):
        """"""
        element_cache_key = self.mesh.elements.___Pr_EWC_cache_key___
        if self.k in (0, 2):
            return element_cache_key
        else: #this is because we have Interface_Dofs_Topology for the dofs of 1-form.
            if  element_cache_key == 'no_cache':
                return 'no_cache'
            else:
                return self.___Pr_EWC_high_similarity_key_generator___

    def ___Pr_EWC_high_similarity_key_generator___(self, i):
        """"""
        element_key = self.mesh.elements.___Pr_EWC_high_similarity_key_generator___(i)
        # noinspection PyUnresolvedReferences
        tT = self.IDT.transition_types # only for 1-form. For 0- and 2-forms, use the key generator of mesh elements.
        if i in tT:
            return element_key + '-' + tT[i]
        else:
            return element_key

    @property
    def orientation(self):
        return self._orientation_

    @property
    def k(self):
        """I am a `k`-form."""
        return self._k_

    @property
    def p(self):
        return self.space.p

    @property
    def coboundary(self):
        return self._coboundary_

    @property
    def numbering(self):
        return self._numbering_

    @property
    def IS(self):
        return self._IS_

    @property
    def BC(self):
        if self._BC_ is None:
            assert self.k != 2, f"I am a volume form."
            self._BC_ = miUsTriangle_SF_BC(self)
        return self._BC_

    def __sub__(self, other):
        """"""
        assert other.__class__.__name__ == self.__class__.__name__
        assert self.mesh == other.mesh
        assert self.space == other.space
        shadow = self.shadow
        COCHAIN = dict()
        for i in self.mesh.elements:
            COCHAIN[i] = self.cochain.local[i] - other.cochain.local[i]
        shadow.cochain.local = COCHAIN
        return shadow

    def __add__(self, other):
        """"""
        assert other.__class__.__name__ == self.__class__.__name__
        assert self.mesh == other.mesh
        assert self.space == other.space
        shadow = self.shadow
        COCHAIN = dict()
        for i in self.mesh.elements:
            COCHAIN[i] = self.cochain.local[i] + other.cochain.local[i]
        shadow.cochain.local = COCHAIN
        return shadow



if __name__ == "__main__":
    # mpiexec -n 4 python
    pass
