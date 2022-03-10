
from inheriting.CSCG.forms.standard.cochain import CSCG_Standard_Form_Cochain_BASE
from _3dCSCG.forms.standard.base.cochain.partial import _3dCSCG_Standard_Form_Cochain_Partial



class _3dCSCG_Standard_Form_Cochain(CSCG_Standard_Form_Cochain_BASE):
    """The cochain must be full. So it must represent all dofs. For partial dofs, we can access them
    through `cochain.partial` property.
    """
    def __init__(self, sf):
        self._partial_ = None
        super().__init__(sf)

    def local_(self, axis):
        """
        The local cochain along a particular axis.

        :param str axis: The local cochain along which axis? ``x``, ``y`` or ``z``.
        :return: The local cochain dict.
        :rtype: Dict[int, numpy.ndarray]
        """
        assert self._sf_.k not in (0, 3), \
            " <Cochain> : %r is a scalar form; has no idea of axes, use cochain.globe." % self._sf_
        numOfBasisComponents = self._sf_.num.basis_components
        localAlongAxis = dict()
        for i in self._sf_.mesh.elements:
            if axis == 'x':
                localAlongAxis[i] = self.local[i][:numOfBasisComponents[0]]
            elif axis == 'y':
                localAlongAxis[i] = self.local[i][numOfBasisComponents[0]:
                                                  numOfBasisComponents[0]+numOfBasisComponents[1]]
            elif axis == 'z':
                localAlongAxis[i] = self.local[i][-numOfBasisComponents[2]:]
            else:
                raise Exception()
        return localAlongAxis


    def ___PRIVATE_local_on_axis___(self, axis, i):
        """
        find the local cochain along a particular axis in a particular mesh_element

        :param str axis: The local cochain along which axis? ``x``, ``y`` or ``z``.
        :param int i: in this mesh element
        :return: The local cochain dict.
        :rtype: Dict[int, numpy.ndarray]
        """
        numOfBasisComponents = self._sf_.num.basis_components
        if axis == 'x':
            localAlongAxis = self.local[i][:numOfBasisComponents[0]]
        elif axis == 'y':
            localAlongAxis = self.local[i][numOfBasisComponents[0]:
                        numOfBasisComponents[0]+numOfBasisComponents[1]]
        elif axis == 'z':
            localAlongAxis = self.local[i][-numOfBasisComponents[2]:]
        else:
            raise Exception()
        return localAlongAxis


    @property
    def partial(self):
        """To access partial of the cochain."""
        if self._partial_ is None:
            self._partial_ = _3dCSCG_Standard_Form_Cochain_Partial(self)
        return self._partial_

