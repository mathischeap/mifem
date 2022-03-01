
import sys
if './' not in sys.path: sys.path.append('./')
from root.config import *
from _2dCSCG.APP.exact_solution.main import ExactSolution
from importlib import import_module
from screws.frozen import FrozenOnly
from _2dCSCG.mesh.domain.inputs.allocator import DomainInputFinder
from _2dCSCG.mesh.domain.main import _2dCSCG_Domain
from _2dCSCG.mesh.main import _2dCSCG_Mesh
from copy import deepcopy
from screws.miscellaneous import MyTimer



class MeshGenerator(FrozenOnly):
    def __init__(self, ID, **kwargs):
        """Remember, **kwargs are parameters to customize the domain.
        The rule is: they can not change the topology of the regions!
        """
        cOmm.barrier()  # for safety reason
        di = DomainInputFinder(ID)(**kwargs)
        self._domain_ = _2dCSCG_Domain(di)
        self._ID_ = ID
        self._kwargs_ = kwargs
        self._freeze_self_()

    @classmethod
    def ___coded_meshes___(cls):
        return DomainInputFinder.___defined_DI___()

    def __call__(self, element_layout, EDM=None, show_info=False):
        if show_info and rAnk == mAster_rank:
            print(f"---[2dCSCG]-[MESH]-{MyTimer.current_time()}-----")
            print(f"   <domain ID>: {self._ID_}")
            str_dp = str(self._kwargs_)
            if len(str_dp) > 40: str_dp = str_dp[:40] + '...'
            print( "   <domain_parameters>: {}".format(str_dp))
            print(f"   <EDM>: {EDM}", flush=True)

        cOmm.barrier()  # for safety reason
        mesh = _2dCSCG_Mesh(self._domain_, element_layout, EDM=EDM)
        mp = dict()
        mp['type'] = '_2dCSCG_Mesh'
        mp['ID'] = self._ID_
        mp['domain_parameters'] = self._kwargs_
        mp['element_layout'] = element_layout
        mp['EDM'] = EDM
        mesh.___define_parameters___ = mp
        dp = dict()
        dp['ID'] = self._ID_
        for key in self._kwargs_:
            dp[key] = self._kwargs_[key]
        mesh.domain.___define_parameters___ = dp
        cOmm.barrier()  # for safety reason

        return mesh



class SpaceInvoker(FrozenOnly):
    def __init__(self, ID):
        cOmm.barrier()  # for safety reason
        assert ID in self.___defined_spaces___(), \
            " <SpaceInvoker> : space <{}> is not coded yet.".format(ID)
        self._ID_ = ID
        cls_name = self.___defined_spaces___()[ID]
        cls_path = self.___space_path___() + ID
        self._space_ = getattr(import_module(cls_path), cls_name)
        self._freeze_self_()

    def __call__(self, inputs, ndim=None, show_info=False):
        if show_info and rAnk == mAster_rank:
            print(f"---[2dCSCG]-[SPACE]-{MyTimer.current_time()}-----")
            print(f"   <space ID>:{self._ID_}")
            print(f"   <space inputs>:{inputs}")
            print(f"   <space ndim>:{ndim}", flush=True)

        cOmm.barrier()  # for safety reason
        if ndim is not None: assert ndim == 2
        S = self._space_(inputs, ndim)
        sp = dict()
        sp['type'] = '_2dCSCG_Space'
        sp['ID'] = self._ID_
        sp['inputs'] = inputs
        sp['ndim'] = ndim
        S.___define_parameters___ = sp
        cOmm.barrier()  # for safety reason
        return S

    @classmethod
    def ___defined_spaces___(cls):
        return {'polynomials': "_2dCSCG_PolynomialSpace"}

    @classmethod
    def ___space_path___(cls):
        return "_2dCSCG.spaces."




class FormCaller(FrozenOnly):
    """Generate a form."""
    def __init__(self, mesh, space):
        cOmm.barrier()  # for safety reason
        self._mesh_ = mesh
        self._space_ = space
        self._freeze_self_()

    def __call__(self, ID, *args, **kwargs):
        cOmm.barrier()  # for safety reason
        cls_path, cls_name = self.___coded_forms___()[ID].split(' : ')
        cls_body = getattr(import_module(cls_path), cls_name)
        if ID in ('scalar', 'vector'):
            FM = cls_body(self._mesh_, *args, **kwargs)
            # continuous forms has no parameters, not a FrozenClass. SO we CAN NOT save them.

        else:
            FM = cls_body(self._mesh_, self._space_, **kwargs)
            fp = dict()
            fp['type'] = '_2dCSCG_Form'
            fp['ID'] = ID
            fp['mesh_parameters'] = deepcopy(self._mesh_.standard_properties.parameters)
            fp['space_parameters'] = deepcopy(self._space_.standard_properties.parameters)
            fp['kwargs'] = kwargs
            FM.___define_parameters___ = fp
        cOmm.barrier()  # for safety reason
        return FM

    @classmethod
    def ___coded_forms___(cls):
        form_path = '_2dCSCG.forms.'
        return {'0-f-i': form_path + "standard._0_form.inner : _2dCSCG_0Form_Inner", # d on inner 0-form is grad
                '1-f-i': form_path + "standard._1_form.inner : _2dCSCG_1Form_Inner", # d on inner 1-form is rot (or curl)
                '2-f-i': form_path + "standard._2_form.inner : _2dCSCG_2Form_Inner",

                '0-f-o': form_path + "standard._0_form.outer : _2dCSCG_0Form_Outer", # d on outer 0-form is curl (or rot)
                '1-f-o': form_path + "standard._1_form.outer : _2dCSCG_1Form_Outer", # d on outer 1-form is div
                '2-f-o': form_path + "standard._2_form.outer : _2dCSCG_2Form_Outer",

                'scalar': "_2dCSCG.fields.scalar : _2dCSCG_ScalarField",
                'vector': "_2dCSCG.fields.vector : _2dCSCG_VectorField",

                '1-t-o': form_path + "trace._1_trace.outer : _2dCSCG_1Trace_Outer",
                }




class ExactSolutionSelector(FrozenOnly):
    """We select an exact solution object with this class."""
    def __init__(self, mesh):
        cOmm.barrier()  # for safety reason
        self._mesh_ = mesh
        self._freeze_self_()

    def __call__(self, ID, show_info=False, **kwargs):
        if show_info and rAnk == mAster_rank:
            print(f"---[2dCSCG]-[Exact Solution]-{MyTimer.current_time()}-----")
            print(f"   <ES ID>: {ID}")
            print(f"   <ES kwargs>: {kwargs}", flush=True)

        cOmm.barrier()  # for safety reason
        assert ID in self.___coded_exact_solution___(), f"Exact solution ID={ID} not found."
        pathAndName = self.___coded_exact_solution___()[ID].split('.')

        classPath = pathAndName[:-1]
        classPath = self.___exact_solution_path___() + '.'.join(classPath)

        className = pathAndName[-1]

        ES =  ExactSolution(self._mesh_)
        status = getattr(import_module(classPath), className)(ES, **kwargs)
        ES.___set_status___(status)
        esp = dict()
        esp['type'] = '_2dCSCG_ExactSolution'
        esp['ID'] = ID
        esp['mesh_parameters'] = deepcopy(self._mesh_.standard_properties.parameters)
        esp['kwargs'] = kwargs
        ES.___define_parameters___ = esp
        cOmm.barrier()  # for safety reason
        return ES

    @classmethod
    def ___coded_exact_solution___(cls):
        return {'sL:sincos1': 'scalar_Laplace.SinCos.SinCos1',}

    @classmethod
    def ___exact_solution_path___(cls):
        return '_2dCSCG.APP.exact_solution.status.'





if __name__ == "__main__":
    # mpiexec python _2dCSCG\main.py

    # mesh = MeshGenerator('cic',)([14,14])
    # mesh = MeshGenerator('crazy',)([14,14])
    mesh = MeshGenerator('cic')([3,3], show_info=True, EDM='chaotic')


    # print(mesh.___PRIVATE_element_division_and_numbering_quality___())
    mesh.visualize.matplot.element_division()
    # mesh.visualize.matplot()

    # print(rAnk, mesh._element_indices_)

    # space = SpaceInvoker('polynomials')([([-1,0.1,1],), ('Lobatto',1)], show_info=True)
    #
    # # from mifem import save, read
    #
    # # mesh.visualize()
    #
    # es = ExactSolutionSelector(mesh)("sL:sincos1", show_info=True)