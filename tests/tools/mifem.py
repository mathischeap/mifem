# -*- coding: utf-8 -*-
"""Here we do tests for the mifem module."""
import sys
if './' not in sys.path:
    sys.path.append('./')
import os
from root.config.main import *

from objects.CSCG._2d.master import MeshGenerator as _2dCSCG_MeshGenerator
from objects.CSCG._2d.master import SpaceInvoker as _2dCSCG_SpaceInvoker
from objects.CSCG._2d.master import FormCaller as _2dCSCG_FormCaller
from objects.CSCG._2d.master import ExactSolutionSelector as _2dCSCG_ExactSolutionSelector

from objects.CSCG._3d.master import MeshGenerator as _3dCSCG_MeshGenerator
from objects.CSCG._3d.master import SpaceInvoker as _3dCSCG_SpaceInvoker
from objects.CSCG._3d.master import ExactSolutionSelector as _3dCSCG_ExactSolutionSelector
from objects.CSCG._3d.master import FormCaller as _3dCSCG_FormCaller
from root.save import save, read

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
    if RANK == MASTER_RANK:
        print("--- [test_mifem_NO1_2dCSCG_save_read] ...... ", flush=True)

    mesh = _2dCSCG_MeshGenerator('cic')([2, 2])
    space = _2dCSCG_SpaceInvoker('polynomials')([('Lobatto', 4), ('Lobatto', 3)])
    FC = _2dCSCG_FormCaller(mesh, space)

    scalar = FC('scalar', p_2d)

    f0 = FC('0-f-i', hybrid=True)
    f0.CF = scalar
    f0.CF.current_time = 0
    f0.discretize()

    save(f0, '_2dCSCG_f0i')

    F0 = read('_2dCSCG_f0i')
    for i in mesh.elements:
        np.testing.assert_array_equal(f0.cochain.local[i], F0.cochain.local[i])
    assert f0.mesh == F0.mesh
    assert f0.space == F0.space

    f0i = FC('0-f-i', hybrid=False)
    f1i = FC('1-f-i', hybrid=False)
    f2i = FC('2-f-i', hybrid=True)
    f0o = FC('0-f-o', hybrid=True)
    f1o = FC('1-f-o', hybrid=True)
    f2o = FC('2-f-o', hybrid=False)
    es = _2dCSCG_ExactSolutionSelector(mesh)('sL:sincos1')

    for f0, f1, f2 in zip([f0i, f0o], [f1i, f1o], [f2i, f2o]):
        f0.CF = es.potential
        f0.CF.current_time = 0
        f0.discretize()
        f0_L2_error = f0.error.L()
        save(f0, '_2dCSCG_f0i')
        F0 = read('_2dCSCG_f0i')
        F0.CF = es.potential
        F0.CF.current_time = 0
        np.testing.assert_almost_equal(F0.error.L(), f0_L2_error)

        f1.CF = es.velocity
        f1.CF.current_time = 0
        f1.discretize()
        f1_L2_error = f1.error.L()
        save(f1, '_2dCSCG_f1i')
        F1 = read('_2dCSCG_f1i')
        F1.CF = es.velocity
        F1.CF.current_time = 0
        np.testing.assert_almost_equal(F1.error.L(), f1_L2_error)

        f2.CF = es.source_term
        f2.CF.current_time = 0
        f2.discretize()
        f2_L2_error = f2.error.L()
        save(f2, '_2dCSCG_f2i')
        F2 = read('_2dCSCG_f2i')
        F2.CF = es.source_term
        F2.CF.current_time = 0
        np.testing.assert_almost_equal(F2.error.L(), f2_L2_error)

        save([f0, f1, f2], 'f0f1f2i')
        F0, F1, F2 = read('f0f1f2i')
        F0.CF = es.potential
        F0.CF.current_time = 0
        F1.CF = es.velocity
        F1.CF.current_time = 0
        F2.CF = es.source_term
        F2.CF.current_time = 0
        np.testing.assert_almost_equal(F0.error.L(), f0_L2_error)
        np.testing.assert_almost_equal(F1.error.L(), f1_L2_error)
        np.testing.assert_almost_equal(F2.error.L(), f2_L2_error)

    if RANK == MASTER_RANK:
        os.remove('_2dCSCG_f0i.mi')
        os.remove('_2dCSCG_f1i.mi')
        os.remove('_2dCSCG_f2i.mi')
        os.remove('f0f1f2i.mi')
    return 1


def test_mifem_NO2_3dCSCG_save_read():
    if RANK == MASTER_RANK:
        print("--- [test_mifem_NO2_3dCSCG_save_read] ...... ", flush=True)

    # ... 3d .......................................

    mesh = _3dCSCG_MeshGenerator('crazy', c=0.0)([2, 1, 1])
    meSH = _3dCSCG_MeshGenerator('crazy', c=0.1)([2, 2, 1])
    save(mesh, '_3dCSCG_mesh')
    MESH = read('_3dCSCG_mesh')
    assert MESH == mesh
    COMM.barrier()
    if RANK == MASTER_RANK:
        os.remove('_3dCSCG_mesh.mi')

    space = _3dCSCG_SpaceInvoker('polynomials')([('Lobatto', 2), ('Lobatto', 1), ('Lobatto', 3)])
    save(space, '_3dCSCG_space')
    SPACE = read('_3dCSCG_space')
    assert SPACE == space
    COMM.barrier()
    if RANK == MASTER_RANK:
        os.remove('_3dCSCG_space.mi')

    es = _3dCSCG_ExactSolutionSelector(MESH)('icpsNS:sincosRD')

    COMM.barrier()

    FC = _3dCSCG_FormCaller(MESH, SPACE)
    vector = FC('vector', (u, v, w))
    f1 = FC('1-f')
    f1.CF = vector
    f1.CF.current_time = 0
    f1.discretize()
    save(f1, '_3dCSCG_1form')
    F1 = read('_3dCSCG_1form')
    COMM.barrier()
    if RANK == MASTER_RANK:
        os.remove('_3dCSCG_1form.mi')
    for i in mesh.elements:
        np.testing.assert_array_equal(f1.cochain.local[i], F1.cochain.local[i])
    assert f1.mesh == F1.mesh
    assert f1.space == F1.space

    f1.CF.current_time = 10
    f1.discretize()
    save(f1, '_3dCSCG_1form')
    F1 = read('_3dCSCG_1form')
    COMM.barrier()
    if RANK == MASTER_RANK:
        os.remove('_3dCSCG_1form.mi')
    f1_error_L = f1.error.L()
    F1.CF = vector
    assert F1.error.L() == f1_error_L

    t2 = FC('2-t')
    t2.CF = es.pressure
    t2.CF.current_time = 11
    t2.discretize()
    save(t2, '_3dCSCG_2trace')
    T2 = read('_3dCSCG_2trace')

    t1 = FC('1-t')
    t1.CF = es.velocity
    t1.CF.current_time = 13
    t1.discretize()
    save(t1, '_3dCSCG_1trace')
    T1 = read('_3dCSCG_1trace')

    t1_cochain = t1.cochain.local
    T1_cochain = T1.cochain.local

    for i in t1_cochain:
        assert i in T1_cochain, f"trivial check."
        np.testing.assert_array_almost_equal(t1_cochain[i], T1_cochain[i])

    t0 = FC('0-t')
    t0.CF = es.pressure
    t0.CF.current_time = 15
    t0.discretize()
    save(t0, '_3dCSCG_0trace')
    T0 = read('_3dCSCG_0trace')

    t0_cochain = t0.cochain.local
    T0_cochain = T0.cochain.local

    for i in t0_cochain:
        assert i in T0_cochain, f"trivial check."
        np.testing.assert_array_almost_equal(t0_cochain[i], T0_cochain[i])

    COMM.barrier()
    if RANK == MASTER_RANK:
        os.remove('_3dCSCG_2trace.mi')
        os.remove('_3dCSCG_1trace.mi')
        os.remove('_3dCSCG_0trace.mi')

    assert t2.CF.current_time == 15

    for i in mesh.elements:
        np.testing.assert_array_equal(t2.cochain.local[i], T2.cochain.local[i])

    save([mesh, MESH, t2, f1, T2, F1, SPACE, space], 'Some_objects')
    SR_mesh, SR_MESH, SR_t2, SR_f1, SR_T2, SR_F1, SR_SPACE, SR_space = read('Some_objects')

    assert SR_t2.mesh is SR_f1.mesh
    assert SR_t2.mesh is SR_F1.mesh
    assert SR_t2.mesh is SR_MESH
    assert SR_t2.space is SR_F1.space

    COMM.barrier()

    assert SR_mesh == mesh == MESH == SR_MESH
    assert SR_SPACE == SR_space == space == SPACE
    for i in mesh.elements:
        np.testing.assert_array_equal(SR_t2.cochain.local[i], T2.cochain.local[i])
    for i in mesh.elements:
        np.testing.assert_array_equal(SR_T2.cochain.local[i], T2.cochain.local[i])

    if RANK == MASTER_RANK:  # Here we generate a random choice that we do not read which objects.
        num_do_not_read = random.randint(0, 8)
        which_do_not_read = random.sample(range(8), num_do_not_read)
    else:
        which_do_not_read = None
    which_do_not_read = COMM.bcast(which_do_not_read, root=MASTER_RANK)
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

    if RANK == MASTER_RANK:
        os.remove('Some_objects.mi')

    return 1


if __name__ == '__main__':
    # mpiexec -n 5 python __tests__/unittests/mifem.py

    test_mifem_NO2_3dCSCG_save_read()

    #
    # test_mifem_NO3_3dCSCG_OLD_read_V0()
