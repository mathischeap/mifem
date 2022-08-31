

from screws.freeze.base import FrozenOnly
from objects.CSCG.base.forms.base.time_wise.func.do import CSCG_FTWF_DO

# noinspection PyUnresolvedReferences
from objects.CSCG._2d.fields.vector.main import _2dCSCG_VectorField
from objects.CSCG._2d.fields.scalar.main import _2dCSCG_ScalarField
# noinspection PyUnresolvedReferences
from objects.CSCG._3d.fields.vector.main import _3dCSCG_VectorField
from objects.CSCG._3d.fields.scalar.main import _3dCSCG_ScalarField


from screws.functions.time_plus_2d_space.constant import CFGt as CFG_2d
from screws.functions.time_plus_3d_space.constant import CFG  as CFG_3d




class CSCG_Form_TimeWise_Func(FrozenOnly):
    def __init__(self, tw):
        self._tw_ = tw
        self._f_ = tw._f_
        self._ndim_ = self._f_.ndim
        self._DO_ = CSCG_FTWF_DO(self)
        self._body_ = None
        self._ES_ = None
        self._ES_variable_name_ = None
        self._freeze_self_()

    def RESET_cache(self):
        """Clear cache."""
        self._body_ = None
        self._ES_ = None
        self._ES_variable_name_ = None

    def ___DO_set_func_body_as___(self, scalar_vector_or_ES, variable_name=None):
        """
        We use this function such that we can set func from exact_solution object.

        If we want the func can be restored, we have to use this function, and we only store and restore
        func from exact_solution object.

        :param scalar_vector_or_ES: scalar, vector or exact_solution, or other special inputs
        :param variable_name:
        :return:
        """

        # - we check if we have one of following special inputs, if yes, we take care of it ----------1
        if isinstance(scalar_vector_or_ES, (int, float)): #-------------------------------------2
            # if so, we make a constant scalar according to the ndim of the mesh
            f = self._tw_._f_
            mesh = f.mesh
            ndim = mesh.ndim

            if ndim == 2:
                cfg = CFG_2d(scalar_vector_or_ES)()

                sca = _2dCSCG_ScalarField(mesh, cfg,
                    ftype='standard', valid_time=None, name=f'constant{scalar_vector_or_ES}')

            elif ndim == 3:
                cfg = CFG_3d(scalar_vector_or_ES)()

                sca = _3dCSCG_ScalarField(mesh, cfg,
                    ftype='standard', valid_time=None, name=f'constant{scalar_vector_or_ES}')
            else:
                raise NotImplementedError()

            scalar_vector_or_ES = sca

        elif isinstance(scalar_vector_or_ES, (list, tuple)) and \
            all([isinstance(_, (int, float)) for _ in scalar_vector_or_ES]): #----------------2
            # a list or tuple of multiple int or float, we make it a vector of constant components
            L = len(scalar_vector_or_ES)
            f = self._tw_._f_
            mesh = f.mesh
            ndim = mesh.ndim
            assert L == ndim, \
                f"I need a list or tuple of {ndim} number as I am in {ndim}-dimensional space."

            fv = list()

            for i in range(L):
                fv.append(eval(f"CFG_{ndim}d(scalar_vector_or_ES[i])()"))

            scalar_vector_or_ES = eval(f"_{ndim}dCSCG_VectorField(mesh, fv, name='constant-vector')")

        else: #==============================================================================2
            pass

        #----------- regular setting process below -------------------------------------------------1
        if scalar_vector_or_ES.__class__.__name__ in (f'_{self._ndim_}dCSCG_ScalarField',
                                                        f'_{self._ndim_}dCSCG_VectorField',
                                                        f'_{self._ndim_}dCSCG_TensorField'):
            assert variable_name is None
            self.body = scalar_vector_or_ES  # has to pass the checker.
        else:
            # only by setting ES and variable_name, the func can be restored.
            self_body = getattr(scalar_vector_or_ES.status, variable_name)
            self._f_.___PRIVATE_TW_FUNC_body_checker___(self_body)  # has to pass the checker.
            # do not use body.setter because we want to set _ES_ and _ES_variable_name_ attributes.
            self._body_ = self_body
            self._ES_ = scalar_vector_or_ES
            self._ES_variable_name_ = variable_name
        #============================================================================================1

    @property
    def body(self):
        """The function body."""
        return self._body_

    @body.setter
    def body(self, body):
        self._f_.___PRIVATE_TW_FUNC_body_checker___(body)  # has to pass the checker.
        self._body_ = body
        self._ES_ = None
        self._ES_variable_name_ = None

    @property
    def ftype(self):
        """(str) The type of the function body: `self.body`."""
        return self.body.ftype

    @property
    def do(self):
        return self._DO_

