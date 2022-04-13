
import sys
if './' not in sys.path: sys.path.append('./')

from objects.CSCG._3d.ADF.trace.base.main import _3dCSCG_Algebra_DUAL_Trace_Form


class _3dCSCG_T0_ADF(_3dCSCG_Algebra_DUAL_Trace_Form):
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
        super().__init__(3, mesh, space, prime, orientation, name)
        self._k_ = 0
        self.standard_properties.___PRIVATE_add_tag___('3dCSCG_standard_algebra_dual_0trace')
        self.___PRIVATE_reset_cache___()
        self._freeze_self_()




    def ___PRIVATE_reset_cache___(self):
        """"""
        super().___PRIVATE_reset_cache___()






if __name__ == "__main__":
    # mpiexec -n 6 python _3dCSCG\ADF\trace\_0_AD_trace.py
    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector

    # mesh = MeshGenerator('bridge_arch_cracked')([8,9,7], EDM='SWV0', show_info=True)
    mesh = MeshGenerator('crazy')([3, 3, 3], EDM=None, show_info=True)
    space = SpaceInvoker('polynomials')([('Lobatto',2), ('Lobatto',1), ('Lobatto',2)], show_info=True)
    # es = ExactSolutionSelector(mesh)('icpsNS:still', show_info=True)
    # ke0 = es.status.kinetic_energy(0)
    # he0 = es.status.helicity(0)


    FC = FormCaller(mesh, space)
    df0 = FC('0-adt', numbering_parameters={'scheme_name': 'Naive', })


    es = ExactSolutionSelector(mesh)('icpsNS:sincosRD')

    df0.prime.TW.func.do.set_func_body_as(es, 'pressure')
    df0.prime.TW.current_time = 0
    df0.prime.TW.do.push_all_to_instant()
    df0.prime.do.discretize()

    # print(df0.cochain.local)

    c0 = df0.cochain.local

    df0.cochain.local = c0