"""Here we do tests for the mifem module."""

import sys
if './' not in sys.path: sys.path.append('./')
import os
from root.config import *

from _2dCSCG.main import MeshGenerator as _2dCSCG_MeshGenerator
from _2dCSCG.main import SpaceInvoker as _2dCSCG_SpaceInvoker
from _2dCSCG.main import FormCaller as _2dCSCG_FormCaller
from _2dCSCG.main import ExactSolutionSelector as _2dCSCG_ExactSolutionSelector

from _3dCSCG.main import MeshGenerator as _3dCSCG_MeshGenerator
from _3dCSCG.main import SpaceInvoker as _3dCSCG_SpaceInvoker
from _3dCSCG.main import ExactSolutionSelector as _3dCSCG_ExactSolutionSelector
from _3dCSCG.main import FormCaller as _3dCSCG_FormCaller
from root.mifem import save, read

import random

def u(t, x, y, z):
    return np.sin(2 * np.pi * x) * np.cos(2 * np.pi * y) * np.cos(2 * np.pi * z) + t / 2

def v(t, x, y, z):
    return np.cos(2 * np.pi * x) * np.sin(2 * np.pi * y) * np.cos(2 * np.pi * z) + t / 2

def w(t, x, y, z):
    return np.cos(2 * np.pi * x) * np.cos(2 * np.pi * y) * np.sin(2 * np.pi * z) + t / 2

def p_2d(t, x, y): return np.cos(2*np.pi*x) * np.cos(2*np.pi*y) + t/2


def test_mifem_NO1_2dCSCG_save_read():
    """
    Unittests for the mesh.
    """
    if rAnk == mAster_rank:
        print("--- [test_mifem_NO1_2dCSCG_save_read] ...... ", flush=True)

    mesh = _2dCSCG_MeshGenerator('cic')([2,2])
    space = _2dCSCG_SpaceInvoker('polynomials')([('Lobatto',4), ('Lobatto',3)])
    FC = _2dCSCG_FormCaller(mesh, space)

    scalar = FC('scalar', p_2d)

    f0 = FC('0-f-i', is_hybrid=True)
    f0.TW.func.___DO_set_func_body_as___(scalar)
    f0.TW.___DO_push_all_to_instant___(0)
    f0.discretize()

    save(f0, '_2dCSCG_f0i')


    F0 = read('_2dCSCG_f0i')
    for i in mesh.elements:
        np.testing.assert_array_equal(f0.cochain.local[i], F0.cochain.local[i])
    assert f0.mesh == F0.mesh
    assert f0.space == F0.space

    f0i = FC('0-f-i', is_hybrid=False)
    f1i = FC('1-f-i', is_hybrid=False)
    f2i = FC('2-f-i', is_hybrid=True)
    f0o = FC('0-f-o', is_hybrid=True)
    f1o = FC('1-f-o', is_hybrid=True)
    f2o = FC('2-f-o', is_hybrid=False)
    es = _2dCSCG_ExactSolutionSelector(mesh)('sL:sincos1')

    for f0, f1, f2 in zip([f0i, f0o], [f1i, f1o], [f2i, f2o]):
        f0.TW.func.___DO_set_func_body_as___(es, 'potential')
        f0.TW.___DO_push_all_to_instant___(0)
        f0.discretize()
        f0_L2_error = f0.error.L()
        save(f0, '_2dCSCG_f0i')
        F0 = read('_2dCSCG_f0i')
        np.testing.assert_almost_equal(F0.error.L(), f0_L2_error)

        f1.TW.func.___DO_set_func_body_as___(es, 'velocity')
        f1.TW.___DO_push_all_to_instant___(0)
        f1.discretize()
        f1_L2_error = f1.error.L()
        save(f1, '_2dCSCG_f1i')
        F1 = read('_2dCSCG_f1i')
        np.testing.assert_almost_equal(F1.error.L(), f1_L2_error)

        f2.TW.func.___DO_set_func_body_as___(es, 'source')
        f2.TW.___DO_push_all_to_instant___(0)
        f2.discretize()
        f2_L2_error = f2.error.L()
        save(f2, '_2dCSCG_f2i')
        F2 = read('_2dCSCG_f2i')
        np.testing.assert_almost_equal(F2.error.L(), f2_L2_error)

        save([f0, f1, f2], 'f0f1f2i')
        F0, F1, F2 = read('f0f1f2i')
        np.testing.assert_almost_equal(F0.error.L(), f0_L2_error)
        np.testing.assert_almost_equal(F1.error.L(), f1_L2_error)
        np.testing.assert_almost_equal(F2.error.L(), f2_L2_error)

    if rAnk == mAster_rank:
        os.remove('_2dCSCG_f0i.mi')
        os.remove('_2dCSCG_f1i.mi')
        os.remove('_2dCSCG_f2i.mi')
        os.remove('f0f1f2i.mi')
    return 1


def test_mifem_NO2_3dCSCG_save_read():
    if rAnk == mAster_rank:
        print("--- [test_mifem_NO2_3dCSCG_save_read] ...... ", flush=True)

    # ... 3d .......................................

    mesh = _3dCSCG_MeshGenerator('crazy', c=0.0)([2,1,1])
    meSH = _3dCSCG_MeshGenerator('crazy', c=0.1)([2,2,1])
    save(mesh, '_3dCSCG_mesh')
    MESH = read('_3dCSCG_mesh')
    assert MESH == mesh
    cOmm.barrier()
    if rAnk == mAster_rank: os.remove('_3dCSCG_mesh.mi')

    space = _3dCSCG_SpaceInvoker('polynomials')([('Lobatto',2), ('Lobatto',1), ('Lobatto',3)])
    save(space, '_3dCSCG_space')
    SPACE = read('_3dCSCG_space')
    assert SPACE == space
    cOmm.barrier()
    if rAnk == mAster_rank: os.remove('_3dCSCG_space.mi')

    es = _3dCSCG_ExactSolutionSelector(MESH)('icpsNS:sincosRD')
    save(es, '_3dCSCG_es')
    ES = read('_3dCSCG_es')
    assert es == ES
    cOmm.barrier()
    if rAnk == mAster_rank: os.remove('_3dCSCG_es.mi')

    FC = _3dCSCG_FormCaller(MESH, SPACE)
    vector = FC('vector', (u,v,w))
    f1 = FC('1-f')
    f1.TW.func.body = vector
    f1.TW.___DO_push_all_to_instant___(0)
    f1.discretize()
    save(f1, '_3dCSCG_1form')
    F1 = read('_3dCSCG_1form')
    cOmm.barrier()
    if rAnk == mAster_rank: os.remove('_3dCSCG_1form.mi')
    for i in mesh.elements:
        np.testing.assert_array_equal(f1.cochain.local[i], F1.cochain.local[i])
    assert f1.mesh == F1.mesh
    assert f1.space == F1.space
    f1.TW.func.___DO_set_func_body_as___(ES, 'velocity')
    f1.TW.current_time = 10
    f1.TW.___DO_push_all_to_instant___()
    f1.discretize()
    save(f1, '_3dCSCG_1form')
    F1 = read('_3dCSCG_1form')
    cOmm.barrier()
    if rAnk == mAster_rank: os.remove('_3dCSCG_1form.mi')
    f1_error_L = f1.error.L()
    assert F1.error.L() == f1_error_L

    t2 = FC('2-t')
    t2.TW.func.___DO_set_func_body_as___(es, 'pressure')
    t2.TW.current_time = 11
    t2.TW.___DO_push_all_to_instant___()
    t2.discretize()
    save(t2, '_3dCSCG_2trace')
    T2 = read('_3dCSCG_2trace')

    t1 = FC('1-t')
    t1.TW.func.DO.set_func_body_as(es, 'velocity')
    t1.TW.current_time = 13
    t1.TW.DO.push_all_to_instant()
    t1.discretize()
    save(t1, '_3dCSCG_1trace')
    T1 = read('_3dCSCG_1trace')
    assert T1.TW.current_time == 13

    t1_cochain = t1.cochain.local
    T1_cochain = T1.cochain.local

    for i in t1_cochain:
        assert i in T1_cochain, f"trivial check."
        np.testing.assert_array_almost_equal(t1_cochain[i], T1_cochain[i])

    t0 = FC('0-t')
    t0.TW.func.DO.set_func_body_as(es, 'pressure')
    t0.TW.current_time = 15
    t0.TW.DO.push_all_to_instant()
    t0.discretize()
    save(t0, '_3dCSCG_0trace')
    T0 = read('_3dCSCG_0trace')

    t0_cochain = t0.cochain.local
    T0_cochain = T0.cochain.local
    assert T0.TW.current_time == 15

    for i in t0_cochain:
        assert i in T0_cochain, f"trivial check."
        np.testing.assert_array_almost_equal(t0_cochain[i], T0_cochain[i])

    cOmm.barrier()
    if rAnk == mAster_rank:
        os.remove('_3dCSCG_2trace.mi')
        os.remove('_3dCSCG_1trace.mi')
        os.remove('_3dCSCG_0trace.mi')

    assert t2.TW.current_time == 11
    for i in mesh.elements:
        np.testing.assert_array_equal(t2.cochain.local[i], T2.cochain.local[i])

    save([mesh, MESH, t2, f1, T2, F1, SPACE, space], 'Some_objects')
    SR_mesh, SR_MESH, SR_t2, SR_f1, SR_T2, SR_F1, SR_SPACE, SR_space = read('Some_objects')

    cOmm.barrier()

    assert SR_mesh == mesh == MESH == SR_MESH
    assert SR_SPACE == SR_space == space == SPACE
    assert SR_t2.TW.current_time == 11
    for i in mesh.elements:
        np.testing.assert_array_equal(SR_t2.cochain.local[i], T2.cochain.local[i])
    assert SR_T2.TW.current_time == 11
    for i in mesh.elements:
        np.testing.assert_array_equal(SR_T2.cochain.local[i], T2.cochain.local[i])
    assert SR_F1.error.L() == f1_error_L
    assert SR_f1.error.L() == f1_error_L

    if rAnk == mAster_rank: # Here we generate a random choice that we do not read which objects.
        num_do_not_read = random.randint(0,8)
        which_do_not_read = random.sample(range(8),num_do_not_read)
    else:
        which_do_not_read = None
    which_do_not_read = cOmm.bcast(which_do_not_read, root=mAster_rank)
    read_individuals = [1 for _ in range(8)]
    for _ in which_do_not_read:
        read_individuals[_] = 0

    READS = read('Some_objects', read_individuals=read_individuals)
    for _ in range(8):
        if _ in which_do_not_read:
            assert READS[_] is None
        else:
            assert READS[_] is not None

    save([mesh, meSH, T1, T0], 'Some_objects')
    read('Some_objects')


    if rAnk == mAster_rank: os.remove('Some_objects.mi')

    return 1


def test_mifem_NO3_3dCSCG_OLD_read_V0():
    if rAnk == mAster_rank:
        print("ooo [test_mifem_NO3_3dCSCG_OLD_read_V0] ...... ", flush=True)

        a = random.random()
    else:
        a = None

    a = cOmm.bcast(a, root=mAster_rank)

    if a > 0.5: # 50% chance to do this test

        u1, u2, w1, w2 = read('TESTS/__unittest_files__/___test_mifem_NO3_read_4_3d_standard_forms___.mi')

        assert u1.mesh._EDM_ is None, "This file is saved at old age, we can only use 'debug' (sequence) EDM."
        assert u2.mesh._EDM_ is None, "This file is saved at old age, we can only use 'debug' (sequence) EDM."
        assert w1.mesh._EDM_ is None, "This file is saved at old age, we can only use 'debug' (sequence) EDM."
        assert w2.mesh._EDM_ is None, "This file is saved at old age, we can only use 'debug' (sequence) EDM."

        if 100 in u1.cochain.local:
            np.testing.assert_array_almost_equal(u1.cochain.local[100],
                         np.array([-7.03038053e-03, -1.00287018e-02, -5.06046544e-03, -7.71643456e-03,
                                   -1.10258838e-02, -5.54211193e-03, -8.72384196e-03, -1.24837939e-02,
                                   -6.25680093e-03, -9.26649387e-03, -1.32843585e-02, -6.67423742e-03,
                                   -7.36697724e-03, -1.06622801e-02, -5.48769521e-03, -8.00219310e-03,
                                   -1.15737091e-02, -5.93213096e-03, -8.95313508e-03, -1.29423681e-02,
                                   -6.59837640e-03, -9.42385645e-03, -1.36118492e-02, -6.93092839e-03,
                                   -8.29207669e-03, -1.22309493e-02, -6.54630560e-03, -8.80267453e-03,
                                   -1.29861633e-02, -6.90160577e-03, -9.43819306e-03, -1.38905988e-02,
                                   -7.29049063e-03, -9.68452517e-03, -1.42245113e-02, -7.41898323e-03,
                                   -8.50618109e-03, -1.25108394e-02, -6.67801526e-03, -8.94248215e-03,
                                   -1.31780867e-02, -7.01584611e-03, -9.45490251e-03, -1.40044663e-02,
                                   -7.43440979e-03, -9.59073419e-03, -1.42945800e-02, -7.59313816e-03,
                                   -4.99085604e-03, -5.96086914e-03, -7.36415935e-03, -8.07106362e-03,
                                   -6.90053217e-03, -8.37554304e-03, -1.05139725e-02, -1.15657646e-02,
                                   -3.38809299e-03, -4.19343120e-03, -5.34204722e-03, -5.88793471e-03,
                                   -4.30192615e-03, -5.29277910e-03, -6.79270414e-03, -7.58714474e-03,
                                   -5.76659734e-03, -7.30939258e-03, -9.59871936e-03, -1.07767773e-02,
                                   -2.69162298e-03, -3.55288847e-03, -4.80973071e-03, -5.43776202e-03,
                                   -1.60363660e-03, -2.87163721e-03, -4.92374571e-03, -6.07018761e-03,
                                   -1.49636685e-03, -3.55897078e-03, -6.89105464e-03, -8.72351188e-03,
                                   -1.66535830e-04, -1.35806356e-03, -3.26961868e-03, -4.31037175e-03,
                                    5.32313669e-05, -1.48654886e-03, -4.00529948e-03, -5.47122939e-03,
                                    1.07741859e-03, -1.59642641e-03, -5.84204334e-03, -8.18389346e-03,
                                    1.42529377e-03, -1.37460593e-04, -2.68364943e-03, -4.06471992e-03,
                                   -3.09119499e-04, -4.37064025e-04, -6.22112877e-04, -7.12667618e-04,
                                   -3.67315167e-04, -5.29626935e-04, -7.74744580e-04, -9.06210988e-04,
                                   -4.71390044e-04, -6.85274603e-04, -1.01925821e-03, -1.19992979e-03,
                                   -5.49183653e-04, -7.90610703e-04, -1.16509558e-03, -1.36692683e-03,
                                   -2.21097025e-03, -2.86781555e-03, -3.83897482e-03, -4.33093719e-03,
                                   -2.46485348e-03, -3.29592928e-03, -4.59396932e-03, -5.28537841e-03,
                                   -2.93773502e-03, -4.07455042e-03, -5.87790763e-03, -6.84519334e-03,
                                   -3.33132489e-03, -4.65034405e-03, -6.69869061e-03, -7.76910486e-03,
                                   -2.40499456e-03, -3.00362740e-03, -3.88992964e-03, -4.33493891e-03,
                                   -2.60289266e-03, -3.41364802e-03, -4.67088729e-03, -5.34554604e-03,
                                   -3.01214173e-03, -4.20246643e-03, -6.08812991e-03, -7.08795741e-03,
                                   -3.36364796e-03, -4.76916173e-03, -6.96667280e-03, -8.10794373e-03]))

        if 8231 in u1.cochain.local:
            np.testing.assert_array_almost_equal(u1.cochain.local[8231],
                         np.array([ 0.00161173,  0.0015913 ,  0.00026379,  0.00232723,  0.00237119,
                                    0.00040466,  0.00317001,  0.00323981,  0.00055588,  0.00360895,
                                    0.00366833,  0.00062405,  0.0015909 ,  0.00155759,  0.0002417 ,
                                    0.00244963,  0.00251886,  0.00042647,  0.00347684,  0.00363419,
                                    0.00063621,  0.00404952,  0.00421501,  0.00074007,  0.00163657,
                                    0.00168323,  0.00023606,  0.00266823,  0.00285302,  0.00046805,
                                    0.00415942,  0.00451709,  0.00080503,  0.00510967,  0.00553803,
                                    0.00100496,  0.00156837,  0.00168802,  0.00019455,  0.00284128,
                                    0.00311813,  0.00050491,  0.00473319,  0.00513803,  0.00090182,
                                    0.00591873,  0.00640905,  0.00114184, -0.01682034, -0.02447582,
                                   -0.03870496, -0.04242771, -0.02768754, -0.03870521, -0.06068099,
                                   -0.06681353, -0.01718276, -0.02320265, -0.03605305, -0.03989559,
                                   -0.0183352 , -0.02542688, -0.03888079, -0.04259164, -0.02998579,
                                   -0.03994532, -0.06080963, -0.06706852, -0.01853824, -0.02383926,
                                   -0.03602956, -0.04001457, -0.0206559 , -0.02649833, -0.0386508 ,
                                   -0.04241659, -0.03360532, -0.0415567 , -0.06038437, -0.06687742,
                                   -0.02071209, -0.02477457, -0.03566153, -0.03985933, -0.02185843,
                                   -0.02688948, -0.03821335, -0.04205432, -0.03560779, -0.04230796,
                                   -0.05974553, -0.06639681, -0.02207622, -0.02536198, -0.03529638,
                                   -0.03958391,  0.00960992,  0.00798095,  0.00303234,  0.00109998,
                                    0.00917224,  0.00758803,  0.00240545,  0.00052323,  0.00821137,
                                    0.00664183,  0.00123466, -0.00048798,  0.00748196,  0.00590092,
                                    0.00046828, -0.00111693,  0.01584447,  0.01284517,  0.00453017,
                                    0.0014635 ,  0.0153156 ,  0.01248527,  0.00377849,  0.00066279,
                                    0.01384081,  0.01120635,  0.0022254 , -0.0008388 ,  0.01264887,
                                    0.01007583,  0.0011483 , -0.00179469,  0.00967916,  0.0077515 ,
                                    0.00263339,  0.00076038,  0.00955654,  0.00776459,  0.00240309,
                                    0.00040761,  0.00886476,  0.00723078,  0.00176438, -0.00034259,
                                    0.00818981,  0.00664783,  0.00129576, -0.00080674]))

        if 3541 in w1.cochain.local:
            np.testing.assert_array_almost_equal(w1.cochain.local[3541],
                         np.array([ 6.21729740e-02,  1.39675981e-01,  9.50832865e-02, -1.49694540e-02,
                                   -7.11951865e-02, -3.10430058e-02, -2.17181114e-01, -2.70439458e-01,
                                   -1.12634020e-01, -3.08110479e-01, -4.50925776e-01, -1.95068912e-01,
                                    2.92515331e-02,  1.64334596e-01,  8.92095679e-02, -5.73898560e-02,
                                   -3.82139866e-02, -2.50961116e-02, -3.26733779e-01, -3.34842862e-01,
                                   -1.32815570e-01, -4.78426161e-01, -5.68577160e-01, -2.44229228e-01,
                                    5.83820506e-02,  2.47409155e-01,  9.05632457e-02,  4.63203327e-02,
                                    8.00815120e-02,  4.45857050e-02, -1.54759935e-01, -6.88786361e-02,
                                   -3.63128853e-02, -2.15469611e-01, -2.70349273e-01, -1.21840601e-01,
                                    3.98807177e-02, -1.95758534e-01, -1.17798598e-01,  8.10782254e-02,
                                   -3.05184000e-02, -8.31776624e-02,  7.69936501e-02,  4.32067036e-02,
                                    3.21869542e-02,  2.11909640e-02,  3.04408577e-04,  4.05125932e-02,
                                   -1.93627812e-01, -1.42117344e-01, -1.28560830e-01, -1.41962113e-01,
                                   -2.95481515e-01, -2.02828857e-01, -1.77559192e-01, -1.81300860e-01,
                                   -1.26989998e-01, -1.52450060e-01, -1.12472517e-01, -9.76313199e-02,
                                   -2.79606120e-01, -1.87417790e-01, -1.25216679e-01, -1.46223395e-01,
                                   -5.98357888e-01, -3.08122408e-01, -2.54473980e-01, -2.30087528e-01,
                                   -5.44738931e-01, -3.19760434e-01, -1.58773333e-01, -1.75481184e-01,
                                   -2.33313095e-01, -2.01090700e-01, -1.64128397e-01, -1.03509520e-01,
                                   -6.22256149e-01, -4.87348028e-01, -4.39251182e-01, -3.06368678e-01,
                                   -5.75991518e-01, -3.98928812e-01, -3.29968595e-01, -2.81239578e-01,
                                   -1.86377724e-01, -1.55706575e-01, -4.20532500e-02,  1.73785196e-02,
                                   -6.40311912e-01, -5.41012470e-01, -3.03878573e-01, -1.57780042e-01,
                                   -4.30088262e-01, -4.50252149e-01, -3.36449834e-01, -2.48667729e-01,
                                    8.39701671e-02, -2.15806701e-01, -3.67612469e-01,  7.09172857e-02,
                                    1.77446098e-02, -1.69324151e-01, -2.50591791e-01, -1.55672607e-01,
                                   -7.80590006e-02, -8.98108393e-02, -1.24441479e-01, -2.74467167e-02,
                                   -2.19948973e-02, -4.16391146e-02, -6.53659322e-02, -3.85135968e-02,
                                    1.57059903e-01, -4.34529425e-01, -4.13610812e-01, -1.20252859e-01,
                                   -8.81467499e-02, -3.62359788e-01, -3.14383949e-01, -2.24410619e-01,
                                   -4.71246295e-03, -1.24950969e-01, -2.74502239e-01, -1.23697488e-01,
                                   -3.35898118e-02, -6.46230122e-02, -1.64671554e-01, -1.46691164e-01,
                                   -7.29396857e-02, -1.65248573e-01,  1.96690374e-02, -3.71267152e-02,
                                   -1.28246342e-01, -1.48049119e-01, -9.93674424e-02, -4.04337787e-02,
                                   -3.95004847e-03, -9.14731450e-02, -1.08625324e-01, -1.24030612e-01,
                                    6.38538174e-02,  7.06419192e-03, -1.34597182e-01, -1.65378722e-01]))

        if 10000 in w1.cochain.local:
            np.testing.assert_array_almost_equal(w1.cochain.local[10000],
                         np.array([-0.12596702,  0.19860482, -0.0532737 , -0.04504215, -0.0591599 ,
                                   -0.12489289, -0.02959033, -0.05586807,  0.1721273 ,  0.00301478,
                                    0.27171651,  0.41074136, -0.0087262 , -0.00348469, -0.1488453 ,
                                   -0.10715691, -0.23675754, -0.09056526, -0.10860165,  0.04551507,
                                    0.17823672,  0.08592895,  0.3033371 ,  0.32454307,  0.022979  ,
                                   -0.08019295, -0.19817858, -0.06534454, -0.05366695, -0.10722873,
                                    0.03694876,  0.05273363,  0.03471565,  0.0510874 ,  0.12352754,
                                    0.20982863, -0.16171696, -0.23243556, -0.33666628, -0.00268864,
                                    0.01299754, -0.05835835,  0.06651712,  0.097189  ,  0.03473298,
                                    0.20198891,  0.29707718,  0.13591641, -0.1339157 ,  0.07705628,
                                    0.20776999,  0.29529613,  0.00441711,  0.10318854, -0.05848773,
                                   -0.62539209,  0.20041473, -0.01396015, -0.23433342, -0.64888386,
                                   -0.15432917,  0.18358809,  0.0293172 ,  0.05116285,  0.05849394,
                                    0.16079227, -0.17220552, -0.28410324,  0.03297839, -0.04130337,
                                   -0.22386796, -0.5973229 ,  0.13187216,  0.01810828, -0.09860641,
                                   -0.14407257,  0.06210268,  0.07817403, -0.00261546, -0.06987283,
                                    0.05420572,  0.04468238, -0.01092409, -0.09321515,  0.06010695,
                                    0.09160675, -0.09135385, -0.14287843,  0.14257266,  0.07243626,
                                    0.07280941, -0.06418665,  0.17764644,  0.07321632,  0.13107583,
                                   -0.06780046,  0.16186484,  0.00729175, -0.27806024,  0.8640307 ,
                                    0.12946627, -0.02662532,  0.0642807 ,  0.08832178, -0.03818082,
                                    0.04489191,  0.07072839, -0.10358556,  0.00353238,  0.05047765,
                                   -0.14600096,  0.11796621,  0.00634068, -0.1985581 ,  0.17773189,
                                    1.06440777, -0.1037303 ,  0.09224652,  0.21164501,  0.34377053,
                                    0.16266334,  0.07215398,  0.10267818, -0.01825721, -0.03788998,
                                    0.0181619 ,  0.06129944,  0.09713491, -0.12771028,  0.06165945,
                                    0.19694763,  0.36760406, -0.04029382,  0.01789267,  0.11617341,
                                    0.25839219,  0.00969525,  0.05363198,  0.0700851 ,  0.0500244 ,
                                    0.15286361,  0.15692786,  0.15166651,  0.06406959]))

        if 2311 in u2.cochain.local:
            np.testing.assert_array_almost_equal(u2.cochain.local[2311],
                     np.array([-2.46622389e-03, -2.34929445e-03, -2.07221922e-03, -1.83648796e-03,
                               -3.43031588e-03, -3.31156529e-03, -2.90256494e-03, -2.51224279e-03,
                               -2.02360235e-03, -1.99421878e-03, -1.77235799e-03, -1.55018847e-03,
                               -3.67281394e-03, -3.51467900e-03, -3.10891498e-03, -2.71273889e-03,
                               -5.14651980e-03, -4.90618340e-03, -4.18551661e-03, -3.61068097e-03,
                               -3.20659155e-03, -3.11355388e-03, -2.68549054e-03, -2.27205231e-03,
                               -1.97557120e-03, -1.91473929e-03, -1.67300544e-03, -1.38283863e-03,
                               -2.70825622e-03, -2.62232255e-03, -2.38042360e-03, -2.12514807e-03,
                               -1.81192444e-03, -1.83789693e-03, -1.73525277e-03, -1.50911336e-03,
                                1.45837056e-05,  2.13144978e-05,  1.24210212e-05, -3.93956372e-04,
                               -7.11023704e-04, -4.64709069e-04, -8.40095929e-04, -1.43784768e-03,
                               -8.52025275e-04, -1.07483419e-03, -1.98501186e-03, -1.27641325e-03,
                                1.92820158e-05,  3.58086648e-05,  2.56495119e-05, -6.54730907e-04,
                               -1.24203037e-03, -9.36868226e-04, -1.48075260e-03, -2.96921442e-03,
                               -2.28341052e-03, -2.11618314e-03, -4.33077834e-03, -3.18570041e-03,
                                1.06583790e-05,  2.42029900e-05,  1.49887206e-05, -4.38703652e-04,
                               -1.04479298e-03, -9.56920403e-04, -1.13407498e-03, -2.87794802e-03,
                               -2.51149923e-03, -1.64099343e-03, -3.72851376e-03, -2.89514374e-03,
                                5.54116929e-04,  1.06292124e-03,  7.45256347e-04,  1.27675200e-03,
                                2.08210544e-03,  1.18112459e-03,  9.46128984e-04,  1.37333075e-03,
                                6.50990865e-04,  8.45725347e-04,  1.51818280e-03,  9.86652697e-04,
                                1.60416380e-03,  2.39996305e-03,  1.17813709e-03,  1.15148933e-03,
                                1.69864461e-03,  8.53210685e-04,  1.36161384e-03,  2.39024759e-03,
                                1.55299059e-03,  2.18988815e-03,  3.40656036e-03,  1.94987124e-03,
                                1.69389731e-03,  2.63214794e-03,  1.34206201e-03,  1.75014279e-03,
                                3.21751513e-03,  2.23473430e-03,  2.79934610e-03,  4.99783737e-03,
                                3.24918693e-03,  2.22677802e-03,  3.38006721e-03,  1.49956966e-03]))

        if 6333 in u2.cochain.local:
            np.testing.assert_array_almost_equal(u2.cochain.local[6333],
                     np.array([-9.45081116e-04, -8.66009495e-04, -7.45144733e-04, -6.78364963e-04,
                               -2.32425589e-03, -2.16564650e-03, -1.92703517e-03, -1.79041525e-03,
                               -2.01184236e-03, -1.91566977e-03, -1.76824199e-03, -1.67520754e-03,
                               -1.54261070e-03, -1.43040007e-03, -1.25134486e-03, -1.14674282e-03,
                               -3.94481213e-03, -3.71448484e-03, -3.35283700e-03, -3.13705725e-03,
                               -3.33932485e-03, -3.19648784e-03, -2.96592252e-03, -2.81615734e-03,
                               -9.67933576e-04, -9.07706203e-04, -8.04657022e-04, -7.40146922e-04,
                               -2.52306710e-03, -2.39396003e-03, -2.18208043e-03, -2.04940920e-03,
                               -2.09070343e-03, -2.00784001e-03, -1.87083694e-03, -1.77951407e-03,
                                4.93168571e-04,  8.62752046e-04,  5.76036810e-04,  3.64219185e-04,
                                6.38916069e-04,  4.32044987e-04,  1.29189576e-04,  2.34073094e-04,
                                1.67206843e-04,  3.50346509e-06,  6.93339301e-06,  4.98813307e-06,
                                5.71711532e-04,  9.99353991e-04,  6.69536212e-04,  4.12247432e-04,
                                7.25780436e-04,  4.95948850e-04,  1.40336951e-04,  2.58987803e-04,
                                1.89948870e-04,  3.95227987e-06,  8.59852700e-06,  6.58485395e-06,
                                2.41307697e-04,  4.16159743e-04,  2.80705065e-04,  1.67614460e-04,
                                2.93703313e-04,  2.04068821e-04,  5.41238281e-05,  1.01330140e-04,
                                7.70807022e-05,  1.61299752e-06,  4.06204748e-06,  3.36956166e-06,
                               -9.50873859e-05, -1.27009398e-04, -7.47974702e-05,  1.17312366e-04,
                                1.52463447e-04,  4.76573601e-05,  2.17151415e-04,  2.86027959e-04,
                                1.24098013e-04, -4.52051625e-05, -2.40354400e-05,  2.42627338e-06,
                                1.93740182e-04,  3.18716860e-04,  1.75880478e-04,  2.46672811e-04,
                                3.65750528e-04,  1.93287372e-04,  2.05163426e-06,  7.05070066e-05,
                                7.14209811e-05,  2.35338691e-04,  4.23892047e-04,  2.66115262e-04,
                                2.40230103e-04,  3.85585357e-04,  2.26891961e-04,  1.55316464e-05,
                                8.99288572e-05,  8.35577065e-05,  2.19730409e-04,  4.04406551e-04,
                                2.60444060e-04,  2.09884625e-04,  3.45856757e-04,  2.09287726e-04]),)

        if 7198 in w2.cochain.local:
            np.testing.assert_array_almost_equal(w2.cochain.local[7198],
                     np.array([ 2.35948467e-04,  2.04162728e-04,  1.44494467e-04,  1.13079053e-04,
                                3.29743943e-04,  3.05743952e-04,  2.11084727e-04,  1.60116927e-04,
                               -7.75079034e-06, -1.49122984e-06, -2.74274969e-06, -8.67211020e-07,
                                1.27319645e-03,  1.09685383e-03,  8.09846005e-04,  5.92742852e-04,
                                1.70042557e-03,  1.51821966e-03,  1.12542219e-03,  7.95757134e-04,
                               -1.01523391e-04, -6.35437709e-05, -3.21021620e-05, -3.39332794e-05,
                                1.28838428e-03,  1.11109334e-03,  7.41788444e-04,  5.01890392e-04,
                                1.64426229e-03,  1.46328664e-03,  9.72867439e-04,  6.53466737e-04,
                               -8.54734660e-05, -6.19050515e-05, -6.82722819e-05, -6.39000328e-05,
                                2.60579345e-04,  4.10249141e-04,  2.26912365e-04,  2.57945737e-04,
                                4.04105095e-04,  2.09067589e-04,  1.42865115e-04,  2.66379275e-04,
                                1.40977183e-04,  7.87198473e-05,  1.79350417e-04,  1.04324825e-04,
                                1.42779663e-03,  2.28484437e-03,  1.33665774e-03,  1.42569468e-03,
                                2.12548141e-03,  1.18321782e-03,  8.18153822e-04,  1.19567813e-03,
                                7.29922477e-04,  4.67706416e-04,  7.26456901e-04,  5.22048027e-04,
                                1.48628167e-03,  2.45738704e-03,  1.42553448e-03,  1.45154832e-03,
                                2.29197087e-03,  1.25407372e-03,  7.64440224e-04,  1.31912500e-03,
                                7.69483920e-04,  4.27971944e-04,  8.71362717e-04,  5.59846940e-04,
                               -4.14884604e-03, -7.99643125e-03, -6.05611894e-03, -9.23599080e-03,
                               -1.51478428e-02, -9.44614103e-03, -2.69804433e-03, -4.05648013e-03,
                               -2.23157177e-03, -4.11442669e-03, -7.93061895e-03, -6.00685875e-03,
                               -9.09691018e-03, -1.49154577e-02, -9.32708282e-03, -2.64015862e-03,
                               -3.96819975e-03, -2.19679495e-03, -3.93598212e-03, -7.48424816e-03,
                               -5.63631569e-03, -8.30716341e-03, -1.35928570e-02, -8.54412242e-03,
                               -2.32769084e-03, -3.53042013e-03, -1.98708938e-03, -3.72395784e-03,
                               -6.94952708e-03, -5.22495687e-03, -7.43907967e-03, -1.21295919e-02,
                               -7.74013192e-03, -2.01479097e-03, -3.07629062e-03, -1.78182465e-03]))

        if 13566 in w2.cochain.local:
            np.testing.assert_array_almost_equal(w2.cochain.local[13566],
                     np.array([ 0.00106803,  0.00122424,  0.00148805,  0.00165126,  0.0027885 ,
                                0.00307746,  0.00351129,  0.00372997,  0.00171858,  0.00191307,
                                0.00220647,  0.00235939,  0.00114887,  0.00131075,  0.00157782,
                                0.00173605,  0.00313772,  0.00343583,  0.00386412,  0.00406752,
                                0.00203377,  0.00223542,  0.00253034,  0.00267126,  0.0001944 ,
                                0.00022457,  0.0002755 ,  0.00030977,  0.00058307,  0.00063689,
                                0.00071451,  0.00075664,  0.00039908,  0.00043783,  0.00049165,
                                0.00052036, -0.00111686, -0.00200859, -0.0013606 , -0.00096018,
                               -0.00151275, -0.00092077, -0.00127526, -0.00188569, -0.0010492 ,
                               -0.00151279, -0.00222165, -0.00121523, -0.00122531, -0.0021317 ,
                               -0.00140311, -0.00107417, -0.00167195, -0.00101588, -0.00158623,
                               -0.00239503, -0.0013736 , -0.00191407, -0.00286491, -0.00161718,
                               -0.00023169, -0.00040703, -0.000271  , -0.00022059, -0.00035789,
                               -0.00022267, -0.00031491, -0.00049438, -0.0002989 , -0.00039112,
                               -0.00060575, -0.00035382,  0.00392028,  0.00644391,  0.00394222,
                                0.00075677,  0.00120877,  0.0007553 ,  0.00015579,  0.00030332,
                                0.00020695,  0.0036074 ,  0.00568425,  0.00333918,  0.00078288,
                                0.00114788,  0.00066505,  0.00019883,  0.00034587,  0.00022006,
                                0.00329436,  0.00495744,  0.00279372,  0.00099684,  0.00144266,
                                0.00081937,  0.00032503,  0.00052083,  0.00032273,  0.0032531 ,
                                0.00485736,  0.00271112,  0.00103733,  0.00150153,  0.00085348,
                                0.00036249,  0.00057838,  0.00034894]))

        if rAnk == mAster_rank:
            print("   +  TEST PASSED.", flush=True)
    else:
        if rAnk == mAster_rank:
            print("   -  TEST SKIPPED.", flush=True)

    return 1


def test_mifem_NO4_3dCSCG_read_V1():
    if rAnk == mAster_rank:
        print("%%% [test_mifem_NO4_3dCSCG_read_V1] ...... ", flush=True)

    if rAnk == mAster_rank: # not always test them.
        a = random.random()
        which = random.randint(0,3)
    else:
        a, which = None, None

    a, which = cOmm.bcast([a, which], root=mAster_rank) # 25% chance to test one form.

    if a > 0.75 and which == 0:
        f0 = read('TESTS/__unittest_files__/test_mifem_NO4_read_3dCSCG_F0.mi')
        np.testing.assert_almost_equal(f0.error.L(), 0.010523898546515776)
        if 1 < sIze < 5000:
            assert f0.mesh._EDM_ is not None, "Must use smart distribution."
        del f0

        if rAnk == mAster_rank:
            print("   ~  Read-3dCSCG-0-form test PASSED.", flush=True)

    else:
        if rAnk == mAster_rank:
            print("   ~  Read-3dCSCG-0-form test SKIPPED.", flush=True)

    if a > 0.75 and which == 1:
        f1 = read('TESTS/__unittest_files__/test_mifem_NO4_read_3dCSCG_F1.mi')
        np.testing.assert_almost_equal(f1.error.L(), 0.03896879461780483)
        if 1 < sIze < 5000:
            assert f1.mesh._EDM_ is not None, "Must use smart distribution."
        del f1

        if rAnk == mAster_rank:
            print("   ~  Read-3dCSCG-1-form test PASSED.", flush=True)

    else:
        if rAnk == mAster_rank:
            print("   ~  Read-3dCSCG-1-form test SKIPPED.", flush=True)

    if a > 0.75 and which == 2:
        f2 = read('TESTS/__unittest_files__/test_mifem_NO4_read_3dCSCG_F2.mi')
        np.testing.assert_almost_equal(f2.error.L(), 0.009633445168902648)
        if 1 < sIze < 5000:
            assert f2.mesh._EDM_ is not None, "Must use smart distribution."
        del f2
        if rAnk == mAster_rank:
            print("   ~  Read-3dCSCG-2-form test PASSED.", flush=True)

    else:
        if rAnk == mAster_rank:
            print("   ~  Read-3dCSCG-2-form test SKIPPED.", flush=True)

    if a > 0.75 and which == 3:
        f3 = read('TESTS/__unittest_files__/test_mifem_NO4_read_3dCSCG_F3.mi')
        np.testing.assert_almost_equal(f3.error.L(), 0.008347731754746129)
        if 1 < sIze < 5000:
            assert f3.mesh._EDM_ is not None, "Must use smart distribution."
        del f3
        if rAnk == mAster_rank:
            print("   ~  Read-3dCSCG-3-form test PASSED.", flush=True)

    else:
        if rAnk == mAster_rank:
            print("   ~  Read-3dCSCG-3-form test SKIPPED.", flush=True)

    return 1


def test_mifem_NO5_2dCSCG_read_V1():
    if rAnk == mAster_rank:
        print("%%% [test_mifem_NO5_2dCSCG_read_V1] ...... ", flush=True)

    if rAnk == mAster_rank: # not always test them.
        a = random.random()
        which = random.randint(0,2)
    else:
        a, which = None, None

    a, which = cOmm.bcast([a, which], root=mAster_rank)

    if a > 0.66 and which == 0:
        f0i = read('TESTS/__unittest_files__/test_mifem_NO5_read_2dCSCG_0if.mi')
        np.testing.assert_almost_equal(f0i.error.L(), 1.022371760237681e-06)
        f0o = read('TESTS/__unittest_files__/test_mifem_NO5_read_2dCSCG_0of.mi')
        np.testing.assert_almost_equal(f0o.error.L(), 1.022371760237681e-06)
        if rAnk == mAster_rank:
            print("   ~  Read-2dCSCG-0-form test PASSED.", flush=True)
    else:
        if rAnk == mAster_rank:
            print("   ~  Read-2dCSCG-0-form test SKIPPED.", flush=True)

    if a > 0.66 and which == 1:
        f1i = read('TESTS/__unittest_files__/test_mifem_NO5_read_2dCSCG_1if.mi')
        np.testing.assert_almost_equal(f1i.error.L(), 0.00030179548111759644)
        f1o = read('TESTS/__unittest_files__/test_mifem_NO5_read_2dCSCG_1of.mi')
        np.testing.assert_almost_equal(f1o.error.L(), 0.00016026043344152729)
        if rAnk == mAster_rank:
            print("   ~  Read-2dCSCG-1-form test PASSED.", flush=True)
    else:
        if rAnk == mAster_rank:
            print("   ~  Read-2dCSCG-1-form test SKIPPED.", flush=True)

    if a > 0.66 and which == 2:
        f2i = read('TESTS/__unittest_files__/test_mifem_NO5_read_2dCSCG_2if.mi')
        np.testing.assert_almost_equal(f2i.error.L(), 3.226148146642923e-05)
        f2o = read('TESTS/__unittest_files__/test_mifem_NO5_read_2dCSCG_2of.mi')
        np.testing.assert_almost_equal(f2o.error.L(), 3.226148146642923e-05)
        if rAnk == mAster_rank:
            print("   ~  Read-2dCSCG-2-form test PASSED.", flush=True)
    else:
        if rAnk == mAster_rank:
            print("   ~  Read-2dCSCG-2-form test SKIPPED.", flush=True)

    return 1





if __name__ == '__main__':
    # mpiexec -n 12 python TESTS\unittest_mifem.py

    test_mifem_NO2_3dCSCG_save_read()

    # test_mifem_NO4_3dCSCG_read_V1()
    #
    # test_mifem_NO3_3dCSCG_OLD_read_V0()