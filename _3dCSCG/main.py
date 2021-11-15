# -*- coding: utf-8 -*-
"""
INTRO

@author: Yi Zhang. Created on Thu Jan 23 00:00:12 2020
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands
"""
import sys
if './' not in sys.path: sys.path.append('./')
from _3dCSCG.APP.exact_solutions.main import ExactSolution
from importlib import import_module
from SCREWS.frozen import FrozenOnly
from SCREWS.miscellaneous import MyTimer
from _3dCSCG.mesh.domain.input.finder import DomainInputFinder
from _3dCSCG.mesh.domain.main import _3dCSCG_Domain
from _3dCSCG.mesh.main import _3dCSCG_Mesh
from copy import deepcopy
from root.config import *


class MeshGenerator(FrozenOnly):
    def __init__(self, ID, **kwargs):
        """Remember, **kwargs are parameters to customize the domain.
        The rule is: they can not change the topology of the regions!
        """
        di = DomainInputFinder(ID)(**kwargs)
        self._domain_ = _3dCSCG_Domain(di)
        self._ID_ = ID
        self._kwargs_ = kwargs
        self._freeze_self_()


    def __call__(self, element_layout, EDM=None, show_info=False):
        """

        :param element_layout:
        :param EDM: ``Element Distribution Method``. When
            EDM = 'debug', we use a naive method
            EDM = None, we find a proper one.
            EDM = a string, we use the specific one.
        :param show_info:
        :return:
        """
        if show_info and rAnk == mAster_rank:
            print(f"---[3dCSCG]-[MESH]-{MyTimer.current_time()}-----")
            print(f"   <domain ID>: {self._ID_}")
            str_dp = str(self._kwargs_)
            if len(str_dp) > 40: str_dp = str_dp[:40] + '...'
            print( "   <domain_parameters>: {}".format(str_dp))
            print(f"   <EDM>: {EDM}", flush=True)

        cOmm.barrier()  # for safety reason
        mesh = _3dCSCG_Mesh(self._domain_, element_layout, EDM=EDM)
        mp = dict()
        mp['type'] = '_3dCSCG_Mesh'
        mp['ID'] = self._ID_
        mp['domain_parameters'] = deepcopy(self._kwargs_)
        mp['element_layout'] = element_layout
        mp['EDM'] = EDM
        mesh.___define_parameters___ = mp

        dp = dict()
        dp['ID'] = self._ID_
        for key in self._kwargs_:
            dp[key] = self._kwargs_[key]
        mesh.domain.___define_parameters___ = dp

        if show_info and rAnk == mAster_rank:
            str_element_layout = str(element_layout)
            if len(str_element_layout) < 40:
                print( "   <element_layout input>: {}".format(str_element_layout))
            else:
                print( "   <element_layout input>: {}...".format(str_element_layout[:40]))
            for rn in mesh.elements.layout:
                print(f"   <element_layout>: {rn} {mesh.elements.layout[rn]}")
            print(f"   <total elements>: {mesh.elements.GLOBAL_num}", flush=True)
        cOmm.barrier()  # for safety reason

        return mesh

    @classmethod
    def ___coded_meshes___(cls):
        return DomainInputFinder.___defined_DI___()



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
            print(f"---[3dCSCG]-[SPACE]-{MyTimer.current_time()}-----")
            print(f"   <space ID>: {self._ID_}")
            print(f"   <space inputs>: {inputs}")
            print(f"   <space ndim>: {ndim}", flush=True)

        cOmm.barrier()  # for safety reason
        S = self._space_(inputs, ndim)
        sp = dict()
        sp['type'] = '_3dCSCG_Space'
        sp['ID'] = self._ID_
        sp['inputs'] = inputs
        sp['ndim'] = ndim
        S.___define_parameters___ = sp
        cOmm.barrier()  # for safety reason
        return S

    @classmethod
    def ___defined_spaces___(cls):
        return {'polynomials': "_3dCSCG_PolynomialSpace"}

    @classmethod
    def ___space_path___(cls):
        return "_3dCSCG.space."





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
        if ID in ('scalar', 'vector', 'tensor'):

            FM = cls_body(self._mesh_, *args, **kwargs)
            # We CANNOT (DO NOT) save continuous field instances.

        else:
            fp = dict()
            fp['type'] = '_3dCSCG_Form'
            fp['ID'] = ID
            fp['mesh_parameters'] = deepcopy(self._mesh_.standard_properties.parameters)
            fp['space_parameters'] = deepcopy(self._space_.standard_properties.parameters)
            fp['kwargs'] = deepcopy(kwargs)

            if ID in ('0-adf', '1-adf', '2-adf', '3-adf',  # algebraic dual (standard) forms
                      '0-adt', '1-adt', '2-adt',           # algebraic dual trace forms
                      ):
                # ---------------- make a dual from a prime ------------------------------------
                if len(args) == 1: # if so, we get a prime form, we make dual form from it.
                    assert kwargs == dict(), \
                        f"when make algebraic dual standard form from prime form, " \
                        f"kwargs must be empty, now it is={kwargs}."
                    prime = args[0]
                    pcn = str(prime.__class__)

                    if ID == '0-adf':
                        assert pcn == "<class '_3dCSCG.form.standard._0_form._0Form'>"
                    elif ID == '1-adf':
                        assert pcn == "<class '_3dCSCG.form.standard._1_form._1Form'>"
                    elif ID == '2-adf':
                        assert pcn == "<class '_3dCSCG.form.standard._2_form._2Form'>"
                    elif ID == '3-adf':
                        assert pcn == "<class '_3dCSCG.form.standard._3_form._3Form'>"
                    elif ID == '0-adt':
                        assert pcn == "<class '_3dCSCG.form.trace._0_trace._0Trace'>"
                    elif ID == '1-adt':
                        assert pcn == "<class '_3dCSCG.form.trace._1_trace._1Trace'>"
                    elif ID == '2-adt':
                        assert pcn == "<class '_3dCSCG.form.trace._2_trace._2Trace'>"
                    else:
                        raise Exception(f"ID={ID} do not accept a single prime form instance as input.")

                    assert prime.IS_hybrid, "prime must be hybrid."
                    assert prime.mesh is self._mesh_ and prime.space is self._space_, \
                        "mesh, space do not match." # not just ==, but is!

                    # if we made a dual from a prime, we still need to feed `orientation` and `name` to the dual class.
                    kwargs = dict()
                    kwargs['orientation'] = prime.orientation
                    kwargs['name'] = 'AD_' + prime.standard_properties.name

                #--------- make the prime from the args and kwargs ---------------------------------
                else: # we make the prime form from the args.
                    p_kwargs = dict()
                    for kw in kwargs:
                        if kw == 'name':
                            p_kwargs[kw] = 'PRIME_' +  kwargs[kw]
                        else:
                            p_kwargs[kw] = kwargs[kw]


                    if ID in ('0-adf', '1-adf', '2-adf', '3-adf'):
                        prime_class_ID = ID.split('-')[0] + '-f'
                        prime = self(prime_class_ID, *args, **p_kwargs, is_hybrid=True)


                    elif ID in ('0-adt', '1-adt', '2-adt'): # note that trace forms must be hybrid.
                        prime_class_ID = ID.split('-')[0] + '-t'
                        prime = self(prime_class_ID, *args, **p_kwargs)


                    else:
                        raise Exception()

                # ---- remove `numbering_parameters` for dual class initialize -------------------
                KWARGS = dict()
                for key in kwargs: # we remove `numbering_parameters` form the kwargs cause we don't need it
                                   # form the initialization of the dual form actually.
                    if key == 'numbering_parameters':
                        pass
                    else:
                        KWARGS[key] = kwargs[key]
                # == ABOVE ==================================================================


                FM = cls_body(prime, self._mesh_, self._space_, **KWARGS)

            else: # all other objects like standard forms, trace forms, edge forms, node forms and so on.

                FM = cls_body(self._mesh_, self._space_, **kwargs)

            FM.___define_parameters___ = fp

        cOmm.barrier()  # for safety reason
        return FM



    @classmethod
    def ___coded_forms___(cls):
        form_path = '_3dCSCG.form.'
        algebra_dual_form_path = '_3dCSCG.ADF.'

        return {'3-f': form_path + "standard._3_form : _3Form",
                '2-f': form_path + "standard._2_form : _2Form",
                '1-f': form_path + "standard._1_form : _1Form",
                '0-f': form_path + "standard._0_form : _0Form",

                '0-adf': algebra_dual_form_path + "standard._0_AD_form : _0_Algebra_DUAL_Form",
                '1-adf': algebra_dual_form_path + "standard._1_AD_form : _1_Algebra_DUAL_Form",
                '2-adf': algebra_dual_form_path + "standard._2_AD_form : _2_Algebra_DUAL_Form",
                '3-adf': algebra_dual_form_path + "standard._3_AD_form : _3_Algebra_DUAL_Form",

                '0-adt': algebra_dual_form_path + "trace._0_AD_trace : _0_Algebra_DUAL_Trace",
                '1-adt': algebra_dual_form_path + "trace._1_AD_trace : _1_Algebra_DUAL_Trace",
                '2-adt': algebra_dual_form_path + "trace._2_AD_trace : _2_Algebra_DUAL_Trace",

                'scalar': "_3dCSCG.field.scalar : _3dCSCG_ScalarField",
                'vector': "_3dCSCG.field.vector : _3dCSCG_VectorField",
                'tensor': "_3dCSCG.field.tensor : _3dCSCG_TensorField",

                '0-t': form_path + "trace._0_trace : _0Trace",
                '1-t': form_path + "trace._1_trace : _1Trace",
                '2-t': form_path + "trace._2_trace : _2Trace",

                '0-e': form_path + "edge._0_edge : _0Edge",
                '1-e': form_path + "edge._1_edge : _1Edge",

                '0-n': form_path + "node._0_node : _0Node",
                }


class ExactSolutionSelector(FrozenOnly):
    """We select an exact solution object with this class."""
    def __init__(self, mesh):
        cOmm.barrier()  # for safety reason
        self._mesh_ = mesh
        self._freeze_self_()

    def __call__(self, ID, show_info=False, **kwargs):
        if show_info and rAnk == mAster_rank:
            print(f"---[3dCSCG]-[Exact Solution]-{MyTimer.current_time()}-----")
            print(f"   <ES ID>: {ID}")
            print(f"   <ES kwargs>: {kwargs}", flush=True)

        cOmm.barrier()  # for safety reason
        assert ID in self.___coded_exact_solution___(), f"Exact solution ID={ID} not found."
        pathAndName = self.___coded_exact_solution___()[ID]
        classPath = pathAndName.split('.')[:-1]
        classPath = self.___exact_solution_path___() + '.'.join(classPath)
        className = pathAndName.split('.')[-1]
        ES =  ExactSolution(self._mesh_)
        status = getattr(import_module(classPath), className)(ES, **kwargs)
        ES.___set_status___(status)
        esp = dict()
        esp['type'] = '_3dCSCG_ExactSolution'
        esp['ID'] = ID
        esp['mesh_parameters'] = deepcopy(self._mesh_.standard_properties.parameters)
        esp['kwargs'] = deepcopy(kwargs)
        ES.___define_parameters___ = esp
        cOmm.barrier()  # for safety reason

        return ES

    @classmethod
    def ___coded_exact_solution___(cls):
        return {'icpsNS:TGV1': 'icpsNS.Taylor_Green_vortex.TGV1',
                'icpsNS:sincosRC': 'icpsNS.Sin_Cos.SinCosRebholz_Conservation',

                'icpsNS:sincos_CCBF': 'icpsNS.Sin_Cos.SinCos_Conservation_Conservative_Body_Force',
                'icpsNS:sincos_CCBF1': 'icpsNS.Sin_Cos.SinCos_Conservation_Conservative_Body_Force1',
                'icpsNS:sincos_CCBF_P': 'icpsNS.Sin_Cos.SinCos_Conservation_Conservative_Body_Force_POLYNOMIALS',
                'icpsNS:sincos_CCBF_C': 'icpsNS.Sin_Cos.SinCos_Conservation_Conservative_Body_Force_CONSTANT',

                'icpsNS:sincosRD': 'icpsNS.Sin_Cos.SinCosRebholz_Dissipation',
                'icpsNS:sincosMD': 'icpsNS.Sin_Cos.SinCos_Modified_Dissipation',
                'icpsNS:CUCD1': 'icpsNS.others.Closed_Unit_Cube_Disspation1',
                'icpsNS:CBFx': 'icpsNS.others.Constant_X_direction_body_force',
                'icpsNS:still': 'icpsNS.others.Still',


                'Poisson:sincos1': "Poisson.Sin_Cos.Poisson_SinCos1",

                }

    @classmethod
    def ___exact_solution_path___(cls):
        return '_3dCSCG.APP.exact_solutions.'







if __name__ == "__main__":
    # mpiexec -n 8 python _3dCSCG\main.py

    mesh = MeshGenerator('bridge_arch_cracked')([1,3,[1,2,4,4,4,4,4,2,1]], EDM='chaotic', show_info=True)
    # mesh = MeshGenerator('bridge_arch_cracked')([2,2,2], EDM='SWV0', show_info=True)
    # mesh = MeshGenerator('crazy')([3, 3, 3], EDM='chaotic', show_info=True)

    es = ExactSolutionSelector(mesh)('Poisson:sincos1')

    print(es.standard_properties.name)

    print(es.status.kinetic_energy(0))


    # mesh.visualize()
    # mesh.visualize.matplot.connection()
    # mesh.domain.visualize()
    # mesh.domain.regions.visualize()

    # print(mesh.domain.regions.map)
    # print(mesh.domain.boundaries.names)
