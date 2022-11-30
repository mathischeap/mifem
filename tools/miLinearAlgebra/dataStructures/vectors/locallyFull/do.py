# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly
from tools.elementwiseCache.gathering.chain import GatheringMatrixChaining

class LocallyFullVectorDo(FrozenOnly):
    """"""
    def __init__(self, LFV):
        self._v_ = LFV
        self._freeze_self_()

    def distributed_to(self, *args, chain_method=None):
        """
        Consider this vector represents cochains of some forms in sequence, we can distribute this vector to the forms.

        :param args: the forms.
        :param chain_method: It can be one of:

            1. ``silly`` -- We distribute the values in sequence to *args. This is corresponded to
                the silly chain method in the chain gathering matrix.
            2. ``sequent`` --
            3. A Chain Gathering Matrix --

        :return:
        """
        if chain_method is None: chain_method = 'silly'

        V = self._v_.V  # V now is 1-d array already.
        
        if chain_method == 'silly':

            if args.count(None) == 0:
                indices = 0
                for form in args:
                    GLOBAL_dofs = form.num.global_dofs
                    # Below we make new locally full vector then distribute it.
                    form.cochain.globe = self._v_.__class__(
                        V[indices : indices + GLOBAL_dofs])
                    indices += GLOBAL_dofs

            elif args.count(None) == 1:
                total_shape = self._v_.shape[0]
                dofs = 0
                for form in args:
                    if form is not None:
                        dofs += form.num.global_dofs
                None_dofs = total_shape - dofs
                if None_dofs <= 0:
                    raise Exception()
                indices = 0
                for form in args:
                    if form is None:
                        indices += None_dofs
                    else:
                        GLOBAL_dofs = form.num.global_dofs
                        form.cochain.globe = self._v_.__class__(
                            V[indices : indices + GLOBAL_dofs])
                        indices += GLOBAL_dofs

            else:
                raise Exception("I can only understand zero or one None")

        elif chain_method == 'sequent':
            cgm = GatheringMatrixChaining(*args)(chain_method='sequent')

            assert cgm.num_GMs == len(args)

            GM0 = cgm.GMs[0]
            local_cochain_all = dict()
            for e in GM0:
                dofs = cgm[e]
                local_cochain_all[e] = V[dofs]

            num_basis = list()
            for f in args:
                _ = f.num.basis
                num_basis.append(_)

            current_local_index = 0
            for i, f in enumerate(args):
                local_cochain = dict()
                for e in GM0:
                    local_cochain[e] = \
                        local_cochain_all[e][current_local_index : current_local_index+num_basis[i]]
                current_local_index += num_basis[i]
                f.cochain.local = local_cochain

        elif chain_method.__class__.__name__  == 'Chain_Gathering_Matrix': # we receive a regular Chain_Gathering_Matrix.
            # we will distribute the entries to forms according to this Chain_Gathering_Matrix.
            cgm = chain_method
            assert cgm.num_GMs == len(args)
            
            GM0 = cgm.GMs[0]
            local_cochain_all = dict()
            for e in GM0:
                dofs = cgm[e]
                local_cochain_all[e] = V[dofs]

            num_basis = list()
            for f in args:
                _ = f.num.basis
                num_basis.append(_)

            current_local_index = 0
            for i, f in enumerate(args):
                local_cochain = dict()
                for e in GM0:
                    local_cochain[e] = \
                        local_cochain_all[e][current_local_index : current_local_index+num_basis[i]]
                current_local_index += num_basis[i]
                f.cochain.local = local_cochain


        else:
            raise NotImplementedError(f'distribution method: {chain_method} not coded.')