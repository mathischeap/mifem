
import sys
if './' not in sys.path: sys.path.append('./')

from _3dCSCG.ADF.trace.base import _3dCSCG_Algebra_DUAL_Trace_Form


class _0_Algebra_DUAL_Trace(_3dCSCG_Algebra_DUAL_Trace_Form):
    """
    Dual trace 0-form.

    :param mesh:
    :param space:
    :param orientation:
    :param name:
    """
    def __init__(self, prime, mesh, space, orientation='outer', name=None):
        """"""
        if name is None: name = orientation + '-oriented-0-AD-Trace'
        super().__init__(3, mesh, space, orientation, name)
        self._prime_ = prime
        self._k_ = 0
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_standard_algebra_dual_0trace')
        self.RESET_cache()
        self._freeze_self_()




    def RESET_cache(self):
        """"""
        super().RESET_cache()






if __name__ == "__main__":
    # mpiexec -n 6 python _3dCSCG\ADF\trace\_0_AD_trace.py
    from _3dCSCG.main import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector

    # mesh = MeshGenerator('bridge_arch_cracked')([8,9,7], EDM='SWV0', show_info=True)
    mesh = MeshGenerator('crazy')([3, 3, 3], EDM=None, show_info=True)
    space = SpaceInvoker('polynomials')([('Lobatto',2), ('Lobatto',1), ('Lobatto',2)], show_info=True)
    # es = ExactSolutionSelector(mesh)('icpsNS:still', show_info=True)
    # ke0 = es.status.kinetic_energy(0)
    # he0 = es.status.helicity(0)


    FC = FormCaller(mesh, space)
    df0 = FC('0-adt', numbering_parameters={'scheme_name': 'Naive', })


    es = ExactSolutionSelector(mesh)('icpsNS:sincosRD')

    df0.prime.TW.func.DO.set_func_body_as(es, 'pressure')
    df0.prime.TW.current_time = 0
    df0.prime.TW.DO.push_all_to_instant()
    df0.prime.DO.discretize()

    # print(df0.cochain.local)

    c0 = df0.cochain.local

    df0.cochain.local = c0