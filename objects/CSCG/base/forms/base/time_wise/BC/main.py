from screws.freeze.base import FrozenOnly
from objects.CSCG.base.forms.base.time_wise.BC.do import CSCG_FTWBC_DO


# noinspection PyUnresolvedReferences
from objects.CSCG._2d.fields.vector.main import _2dCSCG_VectorField
from objects.CSCG._2d.fields.scalar.main import _2dCSCG_ScalarField
# noinspection PyUnresolvedReferences
from objects.CSCG._3d.fields.vector.main import _3dCSCG_VectorField
from objects.CSCG._3d.fields.scalar.main import _3dCSCG_ScalarField


from screws.functions.time_plus_2d_space.constant import CFGt as CFG_2d
from screws.functions.time_plus_3d_space.constant import CFG  as CFG_3d

class CSCG_Form_TimeWise_BC(FrozenOnly):
    """Note that the TW.BC of a form will not be saved!!!!."""
    def __init__(self, tw):
        self._tw_ = tw
        self._f_ = tw._f_
        self._ndim_ = self._f_.ndim
        self._DO_ = CSCG_FTWBC_DO(self)
        self._body_ = None
        self._ES_ = None
        self._ES_variable_name_ = None
        self._freeze_self_()

    def RESET_cache(self):
        """Clear cache."""
        self._body_ = None
        self._ES_ = None
        self._ES_variable_name_ = None


    def ___DO_set_BC_body_as___(self, scalar_vector):
        """
        We use this function such that we can set BC.

        :param scalar_vector: scalar, vector
        :return:
        """
        if isinstance(scalar_vector, (int, float)): #------------------------------------------2
            # if so, we make a constant scalar according to the ndim of the mesh
            f = self._tw_._f_
            mesh = f.mesh
            ndim = mesh.ndim

            if ndim == 2:
                cfg = CFG_2d(scalar_vector)

                sca = _2dCSCG_ScalarField(mesh, cfg,
                    ftype='standard', valid_time=None, name=f'constant{scalar_vector}')

            elif ndim == 3:
                cfg = CFG_3d(scalar_vector)

                sca = _3dCSCG_ScalarField(mesh, cfg,
                    ftype='standard', valid_time=None, name=f'constant{scalar_vector}')
            else:
                raise NotImplementedError()

            scalar_vector = sca

        elif isinstance(scalar_vector, (list, tuple)) and \
            all([isinstance(_, (int, float)) for _ in scalar_vector]): #-----------------------2
            # a list or tuple of multiple int or float, we make it a vector of constant components
            L = len(scalar_vector)
            f = self._tw_._f_
            mesh = f.mesh
            ndim = mesh.ndim
            assert L == ndim, \
                f"I need a list or tuple of {ndim} number as I am in {ndim}-dimensional space."

            fv = list()

            for i in range(L):
                fv.append(eval(f"CFG_{ndim}d(scalar_vector[i])"))

            scalar_vector = eval(f"_{ndim}dCSCG_VectorField(mesh, fv, name='constant-vector')")

        else: #================================================================================2
            pass

        return scalar_vector

    @property
    def body(self):
        """The function body."""
        return self._body_

    @body.setter
    def body(self, body):
        body = self.___DO_set_BC_body_as___(body)
        self._f_.___PRIVATE_TW_BC_body_checker___(body) # has to pass the checker.
        self._body_ = body
        self._ES_ = None
        self._ES_variable_name_ = None

    @property
    def ftype(self):
        """(str) The function type."""
        return self.body.ftype

    @property
    def do(self):
        return self._DO_