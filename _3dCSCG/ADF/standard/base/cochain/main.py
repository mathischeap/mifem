

from root.config.main import *
from screws.freeze.main import FrozenOnly
from _3dCSCG.ADF.standard.base.cochain.local import ____3dCSCG_ADSF_Cochain_Local____
from scipy.sparse import csr_matrix
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.main import EWC_ColumnVector




class _3dCSCG_Algebra_DUAL_Standard_Form_Cochain(FrozenOnly):
    """The cochain of algebra dual form is equal to the mass matrix dot the cochain of the prime form.

    dual_cochain = mass_matrix dot prime_cochain.

    """
    def __init__(self, dsf):
        """
        :param dsf: The dual standard form.
        """
        self._dsf_ = dsf
        self._local_ = None  # this is a key property, should not reset it.
        self._freeze_self_()

    @property
    def local(self):
        """We know that the local cochain of the prime form is a dict whose keys are local element numbers and values
        are the local cochains (1-d array). While for algebra dual standard forms, we make a EWC_ColumnVector for it
        since we do not want to save the local cochain of the algebra dual form. We will generate the cochain when
        we call it in real time.

        :return:
        """
        if self._local_ is None:
            self._local_ = ____3dCSCG_ADSF_Cochain_Local____(self)
            # the local cochain will be renewed in real-time if the local cochain of the prime form is renewed.
        return self._local_

    @local.setter
    def local(self, local):
        """"""
        assert len(local) == len(self._dsf_.mesh.elements), f"length of local is wrong!"
        prime_local_cochain = dict()
        iM = self._dsf_.inverse_mass_matrix
        for i in self._dsf_.mesh.elements:
            assert i in local, f"mesh element #{i} is missing in local in Core #{rAnk}."
            prime_local_cochain[i] = iM[i] @ local[i]
        self._dsf_.prime.cochain.local = prime_local_cochain

    @property
    def EWC(self):
        """For dual standard forms, this is actually very similar with the local.

        :return:
        """
        return EWC_ColumnVector(self._dsf_.mesh.elements,
                                self.___PRIVATE_local_cochain_link___,
                                cache_key_generator='no_cache')

    def ___PRIVATE_local_cochain_link___(self, i):
        """"""
        MM = self._dsf_.mass_matrix[i]
        LC = self._dsf_.prime.cochain.local[i]
        DLC = csr_matrix(MM @ LC).T
        return DLC


    @property
    def globe(self):
        raise NotImplementedError()

    @globe.setter
    def globe(self, globe):
        """We have to set the cochain to the local cochain of the prime.

        :param globe:
        :return:
        """
        raise NotImplementedError()

    def __getitem__(self, i):
        """If `i` is an element number in this core, we should be able to return a local cochain of it (if not None)"""
        return self.local[i]

    def __contains__(self, i):
        """Return if an element number (`i`) is in this core."""
        return i in self.local

    def __iter__(self):
        """Go through all mesh element numbers in this core."""
        for i in self.local:
            yield i

    def __len__(self):
        """Actually return how many mesh elements in this core."""
        return len(self.local)