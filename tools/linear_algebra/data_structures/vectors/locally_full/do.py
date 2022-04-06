from screws.freeze.base import FrozenOnly


class LocallyFullVectorDo(FrozenOnly):
    """"""
    def __init__(self, LFV):
        self._v_ = LFV
        self._freeze_self_()

    def distributed_to(self, *args, method='sequence'):
        """
        Consider this vector represents cochains of some forms in sequence, we can distribute this vector to the forms.

        :param args: the forms.
        :param method: It can be one of:

            1. ``sequence`` -- We distribute the values in sequence to *args.

        :return:
        """
        if method == 'sequence':

            V = self._v_.V # V now is 1-d array already.

            if args.count(None) == 0:
                indices = 0
                for form in args:
                    GLOBAL_dofs = form.num.GLOBAL_dofs
                    # Below we make new locally full vector then distribute it.
                    form.cochain.globe = self._v_.__class__(
                        V[indices : indices + GLOBAL_dofs])
                    indices += GLOBAL_dofs

            elif args.count(None) == 1:
                total_shape = self._v_.shape[0]
                dofs = 0
                for form in args:
                    if form is not None:
                        dofs += form.num.GLOBAL_dofs
                None_dofs = total_shape - dofs
                if None_dofs <= 0:
                    raise Exception()
                indices = 0
                for form in args:
                    if form is None:
                        indices += None_dofs
                    else:
                        GLOBAL_dofs = form.num.GLOBAL_dofs
                        form.cochain.globe = self._v_.__class__(
                            V[indices : indices + GLOBAL_dofs])
                        indices += GLOBAL_dofs
            else:
                raise Exception("I can only understand one ")
            # for ag in args:
            #     if ag is None:
            #         pass

        else:
            raise NotImplementedError(f'distribution method: {method} not coded.')