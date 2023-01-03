# -*- coding: utf-8 -*-
from root.config.main import RANK, MASTER_RANK, COMM, np
from components.freeze.main import FrozenOnly
from objects.CSCG._3d.ADF.trace.base.cochain.local import ____3dCSCG_ADTF_Cochain_Local____
from scipy.sparse import csc_matrix


class _3dCSCG_Algebra_DUAL_Trace_Form_Cochain(FrozenOnly):
    """The cochain of algebra dual form is equal to the mass matrix dot the cochain of the prime
    form.

    dual_cochain = mass_matrix dot prime_cochain.

    """
    def __init__(self, dt):
        """
        :param dt: The dual standard form.
        """
        self._dt_ = dt
        self._local_ = None  # this is a key property, should not reset it.
        self._freeze_self_()

    @property
    def local(self):
        """We know that the local cochain of the prime form is a dict whose keys are local element
        numbers and values are the local cochains (1-d array). While for algebra dual standard
        forms, we make a EWC_ColumnVector for it since we do not want to save the local cochain of
        the algebra dual form. We will generate the cochain when we call it in real time.

        :return:
        """
        if self._local_ is None:
            self._local_ = ____3dCSCG_ADTF_Cochain_Local____(self)
            # the local cochain will be renewed automatically if the local cochain of the prime form is renewed.
        return self._local_

    @local.setter
    def local(self, local):
        """"""
        assert len(local) == len(self._dt_.mesh.elements), \
            f"length of local is wrong!"
        prime_local_cochain = dict()
        iM = self._dt_.inverse_mass_matrix
        for i in self._dt_.mesh.elements:
            assert i in local, \
                f"mesh element #{i} is missing in local in Core #{RANK}."
            prime_local_cochain[i] = iM[i] @ local[i]
        self._dt_.prime.cochain.local = prime_local_cochain

    @property
    def globe(self):
        raise NotImplementedError()

    @globe.setter
    def globe(self, globe):
        """We have to set the cochain to the local cochain of the prime.

        :param globe:
        :return:
        """
        if globe.__class__.__name__ == 'DistributedVector':
            assert globe.V.shape == (self._dt_.num.global_dofs, 1), "globe cochain shape wrong."
            # gather vector to master core ...
            if globe.IS_master_dominating:
                # no need to gather
                VV = globe.V.T.toarray()[0]
            else:
                V = globe.V
                V = COMM.gather(V, root=MASTER_RANK)
                if RANK == MASTER_RANK:
                    VV = np.empty((self._dt_.num.global_dofs,))
                    for v in V:
                        indices = v.indices
                        data = v.data
                        VV[indices] = data
            # distribute vector to individual cores ...
            local_range = self._dt_.numbering.gathering.local_range
            local_range = COMM.gather(local_range, root=MASTER_RANK)
            if RANK == MASTER_RANK:
                TO_BE_SENT = list()
                for lr in local_range:
                    if lr == tuple():
                        to_be_sent = None
                    else:
                        # noinspection PyUnboundLocalVariable
                        to_be_sent = csc_matrix(
                            (VV[lr[0]:lr[1]], range(lr[0], lr[1]), [0, lr[1]-lr[0]]),
                            shape=(self._dt_.num.global_dofs, 1))
                    TO_BE_SENT.append(to_be_sent)
            else:
                TO_BE_SENT = None
            TO_BE_SENT = COMM.scatter(TO_BE_SENT, root=MASTER_RANK)
            # distribute to local cochain ...
            local = dict()
            GM = self._dt_.prime.numbering.gathering
            for i in GM:  # go through all local elements
                idx = GM[i].full_vector
                local[i] = TO_BE_SENT[idx].toarray().ravel()
            self.local = local

        elif globe.__class__.__name__ == 'LocallyFullVector':
            V = globe.V  # V already be 1-d array.
            local = dict()
            GM = self._dt_.prime.numbering.gathering
            for i in GM:  # go through all local elements
                idx = GM[i].full_vector
                local[i] = V[idx]
            self.local = local

        else:
            raise Exception(f"Can not set cochain from {globe}.")

    def __getitem__(self, i):
        """If `i` is an element number in this core, we should be able to return a local cochain of
        it (if not None)
        """
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
