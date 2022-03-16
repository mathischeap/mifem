




from inheriting.CSCG.forms.standard.cochain import CSCG_Standard_Form_Cochain_BASE



class _2dCSCG_Standard_Form_Cochain(CSCG_Standard_Form_Cochain_BASE):
    def __init__(self, sf):
        super().__init__(sf)

    def local_(self, axis):
        """
        The local cochain along a particular axis.

        :param str axis: The local cochain along which axis? ``x``, ``y``.
        :return: The local cochain dict.
        :rtype: Dict[int, numpy.ndarray]
        """
        assert self._sf_.k not in (0, 2), \
            " <Cochain> : %r is a scalar form; has no idea of axes, use cochain.globe." % self._sf_
        numOfBasisComponents = self._sf_.num.basis_components
        localAlongAxis = dict()
        for i in self._sf_.mesh.elements:
            if axis == 'x':
                localAlongAxis[i] = self.local[i][:numOfBasisComponents[0]]
            elif axis == 'y':
                localAlongAxis[i] = self.local[i][numOfBasisComponents[0]:]
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
            localAlongAxis = self.local[i][numOfBasisComponents[0]:]
        else:
            raise Exception()
        return localAlongAxis