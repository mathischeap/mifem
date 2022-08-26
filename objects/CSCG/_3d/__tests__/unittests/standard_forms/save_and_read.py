# -*- coding: utf-8 -*-

import os
import sys
if './' not in sys.path: sys.path.append('./')


from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller, ExactSolutionSelector
from screws.miscellaneous.mios import remove
from root.config.main import mAster_rank, rAnk
from root.save import save
from root.read.main import read

def test_standard_forms_save_and_read():
    """"""
    if rAnk == mAster_rank:
        print("SSR [test_standard_forms_save_and_read] ...... ", flush=True)

    ## ---------- mesh and space ------------------------------------------------------------
    mesh = MeshGenerator('bridge_arch_cracked',)([5,4,7], EDM='debug')
    space = SpaceInvoker('polynomials')([('Lobatto', 3), ('Lobatto', 4), ('Lobatto', 2)])
    FC = FormCaller(mesh, space)
    es = ExactSolutionSelector(mesh)('icpsNS:sincosRD')

    ## --------- make the data -------------------------------------------------------------
    # f0 = FC('0-f', is_hybrid=True)
    # f0.TW.func.do.set_func_body_as(es, 'pressure')
    # f0.TW.do.push_all_to_instant(0)
    # f0.discretize()
    # save(f0, '___f0_mi_SFsr___.mi')
    # print(f0.error.L())
    #
    # f1 = FC('1-f', is_hybrid=False)
    # f1.TW.func.do.set_func_body_as(es, 'velocity')
    # f1.TW.do.push_all_to_instant(0)
    # f1.discretize()
    # save(f1, '___f1_mi_SFsr___.mi')
    # print(f1.error.L())
    #
    #
    # f2 = FC('2-f', is_hybrid=True)
    # f2.TW.func.do.set_func_body_as(es, 'velocity')
    # f2.TW.do.push_all_to_instant(0)
    # f2.discretize()
    # save(f2, '___f2_mi_SFsr___.mi')
    # print(f2.error.L())
    #
    # f3 = FC('3-f')
    # f3.TW.func.do.set_func_body_as(es, 'pressure')
    # f3.TW.do.push_all_to_instant(0)
    # f3.discretize()
    # save(f3, '___f3_mi_SFsr___.mi')
    # print(f3.error.L())
    #
    # save([f0, f1, f2, f3], '___f0123_mi_SFsr___.mi')

    absolute_path = os.path.dirname(__file__)
    f = read(absolute_path + '/auxiliaries/___f0_mi_SFsr___.mi')
    # print(f.error.L())
    assert f.error.L() < 0.02
    f = read(absolute_path + '/auxiliaries/___f1_mi_SFsr___.mi')
    assert f.error.L() < 0.02
    f = read(absolute_path + '/auxiliaries/___f2_mi_SFsr___.mi')
    assert f.error.L() < 0.09
    f = read(absolute_path + '/auxiliaries/___f3_mi_SFsr___.mi')
    assert f.error.L() < 0.05

    f0, f1, f2, f3 = read(absolute_path + '/auxiliaries/___f0123_mi_SFsr___.mi')
    assert f0.error.L() < 0.02 and f1.error.L() < 0.02 and f2.error.L() < 0.09 and f3.error.L() < 0.07

    ## --------- make the data -------------------------------------------------------------
    # df1 = FC('1-adf')
    # dt0 = FC('0-adt')
    #
    # df1.prime.TW.func.do.set_func_body_as(es, 'velocity')
    # df1.prime.TW.current_time = 1
    # df1.prime.TW.do.push_all_to_instant()
    # df1.prime.do.discretize()
    #
    # dt0.prime.TW.func.do.set_func_body_as(es, 'pressure')
    # dt0.prime.TW.current_time = 1
    # dt0.prime.TW.do.push_all_to_instant()
    # dt0.prime.do.discretize()
    # save([df1, dt0], '___df1_dt0_mi_SFsr___.mi')

    df1, dt0 = read(absolute_path + '/auxiliaries/___df1_dt0_mi_SFsr___.mi')
    df0 = df1.coboundary(dt0)
    df0.prime.TW.func.do.set_func_body_as(es, 'pressure')
    df0.prime.TW.current_time = 1
    df0.prime.TW.do.push_all_to_instant()
    assert 170 < df0.prime.error.L() < 175

    df2 = FC('2-adf')
    df2.prime.TW.func.do.set_func_body_as(es, 'velocity')
    df2.prime.TW.current_time = 1
    df2.prime.TW.do.push_all_to_instant()
    df2.prime.do.discretize()
    save([df2, df1], '___df2_df1_mi_SFsr___.mi')
    df2, df1 = read('___df2_df1_mi_SFsr___.mi')
    assert df2.mesh is df1.mesh
    remove('___df2_df1_mi_SFsr___.mi')

    return 1


if __name__ == '__main__':
    # mpiexec -n 4 python objects\CSCG\_3d\__tests__\unittests\standard_forms\save_and_read.py
    test_standard_forms_save_and_read()