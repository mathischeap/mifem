# -*- coding: utf-8 -*-
"""
@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path:
    sys.path.append('./')
from root.config.main import *
from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller
from tests.objects.CSCG._3d.randObj.form_caller import random_FormCaller_of_total_load_around
import random



def test_Naive_Numbering_NO1_0form():
    """"""
    if RANK == MASTER_RANK:
        print("--- [test_Naive_Numbering_NO1_0form] ...... ", flush=True)

    mesh = MeshGenerator('crazy_periodic')([2, 2, 2], EDM='debug')
    space = SpaceInvoker('polynomials')([('Lobatto', 2), ('Lobatto', 1), ('Lobatto', 2)])
    FC = FormCaller(mesh, space)
    f0 = FC('0-f', hybrid=False, numbering_parameters='Naive')

    benchmark = np.array(
        [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17],
         [2, 18, 0, 5, 21, 3, 8, 19, 6, 11, 22, 9, 14, 20, 12, 17, 23, 15],
         [3, 4, 5, 0, 1, 2, 9, 10, 11, 6, 7, 8, 15, 16, 17, 12, 13, 14],
         [5, 21, 3, 2, 18, 0, 11, 22, 9, 8, 19, 6, 17, 23, 15, 14, 20, 12],
         [12, 13, 14, 15, 16, 17, 24, 26, 28, 25, 27, 29, 0, 1, 2, 3, 4, 5],
         [14, 20, 12, 17, 23, 15, 28, 30, 24, 29, 31, 25, 2, 18, 0, 5, 21, 3],
         [15, 16, 17, 12, 13, 14, 25, 27, 29, 24, 26, 28, 3, 4, 5, 0, 1, 2],
         [17, 23, 15, 14, 20, 12, 29, 31, 25, 28, 30, 24, 5, 21, 3, 2, 18, 0]]
    )
    for i in f0.numbering.gathering:
        assert np.all(f0.numbering.gathering[i].full_vector == benchmark[i, :])

    return 1


def test_Naive_Numbering_NO2_1form():
    """"""
    if RANK == MASTER_RANK:
        print("--- [test_Naive_Numbering_NO2_1form] ...... ", flush=True)


    mesh = MeshGenerator('crazy_periodic')([2, 2, 2], EDM='debug')
    space = SpaceInvoker('polynomials')([('Lobatto', 2), ('Lobatto', 1), ('Lobatto', 2)])
    FC = FormCaller(mesh, space)
    f1 = FC('1-f', hybrid=False, numbering_parameters='Naive')

    benchmark = np.array(
        [[0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
          17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32],
         [33, 39, 36, 42, 34, 40, 37, 43, 35, 41, 38, 44, 14, 45, 12, 17, 46,
          15, 20, 47, 18, 23, 48, 21, 26, 50, 24, 29, 49, 27, 32, 51, 30],
         [2, 3, 0, 1, 6, 7, 4, 5, 10, 11, 8, 9, 52, 55, 58, 53, 56,
          59, 54, 57, 60, 24, 25, 26, 21, 22, 23, 30, 31, 32, 27, 28, 29],
         [36, 42, 33, 39, 37, 43, 34, 40, 38, 44, 35, 41, 58, 61, 52, 59, 62,
          53, 60, 63, 54, 26, 50, 24, 23, 48, 21, 32, 51, 30, 29, 49, 27],
         [8, 9, 10, 11, 64, 66, 65, 67, 0, 1, 2, 3, 18, 19, 20, 68, 69,
          70, 12, 13, 14, 71, 75, 79, 73, 77, 81, 72, 76, 80, 74, 78, 82],
         [35, 41, 38, 44, 83, 85, 84, 86, 33, 39, 36, 42, 20, 47, 18, 70, 87,
          68, 14, 45, 12, 79, 88, 71, 81, 90, 73, 80, 89, 72, 82, 91, 74],
         [10, 11, 8, 9, 65, 67, 64, 66, 2, 3, 0, 1, 54, 57, 60, 92, 93,
          94, 52, 55, 58, 73, 77, 81, 71, 75, 79, 74, 78, 82, 72, 76, 80],
         [38, 44, 35, 41, 84, 86, 83, 85, 36, 42, 33, 39, 60, 63, 54, 94, 95,
          92, 58, 61, 52, 81, 90, 73, 79, 88, 71, 82, 91, 74, 80, 89, 72]]
    )
    for i in f1.numbering.gathering:
        assert np.all(f1.numbering.gathering[i].full_vector == benchmark[i, :])

    return 1


def test_Naive_Numbering_NO3_2form():
    """"""
    if RANK == MASTER_RANK:
        print("--- [test_Naive_Numbering_NO3_2form] ...... ", flush=True)

    mesh = MeshGenerator('crazy_periodic')([2, 2, 2], EDM='debug')
    space = SpaceInvoker('polynomials')([('Lobatto', 2), ('Lobatto', 1), ('Lobatto', 2)])
    FC = FormCaller(mesh, space)

    f2 = FC('2-f', hybrid=False, numbering_parameters='Naive')
    benchmark = np.array([
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
         17, 18, 19],
        [2, 20, 0, 5, 21, 3, 22, 26, 24, 28, 23, 27, 25, 29, 30, 33, 31,
         34, 32, 35],
        [36, 38, 40, 37, 39, 41, 8, 9, 6, 7, 12, 13, 10, 11, 42, 45, 43,
         46, 44, 47],
        [40, 48, 36, 41, 49, 37, 24, 28, 22, 26, 25, 29, 23, 27, 50, 53, 51,
         54, 52, 55],
        [56, 58, 60, 57, 59, 61, 62, 66, 64, 68, 63, 67, 65, 69, 18, 19, 70,
         71, 14, 15],
        [60, 72, 56, 61, 73, 57, 74, 78, 76, 80, 75, 79, 77, 81, 32, 35, 82,
         83, 30, 33],
        [84, 86, 88, 85, 87, 89, 64, 68, 62, 66, 65, 69, 63, 67, 44, 47, 90,
         91, 42, 45],
        [88, 92, 84, 89, 93, 85, 76, 80, 74, 78, 77, 81, 75, 79, 52, 55, 94,
         95, 50, 53]
    ])
    for i in f2.numbering.gathering:
        assert np.all(f2.numbering.gathering[i].full_vector == benchmark[i, :])

    return 1


# noinspection PyUnboundLocalVariable
def test_Naive_Numbering_NO5_0trace():
    """"""
    if RANK == MASTER_RANK:
        load = random.randint(10, 199)
        print(f"--- [test_Naive_Numbering_NO5_0trace] @ FC-load = {load} ...... ", flush=True)
    else:
        load = None
    load = COMM.bcast(load, root=MASTER_RANK)
    FC = random_FormCaller_of_total_load_around(load)

    t0 = FC('0-t', numbering_parameters={'scheme_name': 'Naive'})

    GM_TEW = t0.numbering.trace_element_wise
    for i in GM_TEW:
        assert i in t0.mesh.trace.elements, "must be the case"
        assert len(set(GM_TEW.keys())) == t0.mesh.trace.elements.num, "must be the case"

    GM_TEW = COMM.gather(GM_TEW, root=MASTER_RANK)
    if RANK == MASTER_RANK:
        END = 0
        for i in range(t0.mesh.trace.elements.global_num):  # go through all (global) trace elements
            for gm_core in GM_TEW:
                if i in gm_core:
                    fv = gm_core[i].full_vector
                    assert fv[0] == END, f"must be the case."
                    NEW_END = fv[-1] + 1

            END = NEW_END

    num_of_dofs_in_this_core = t0.numbering.num_local_dofs
    LOCAL_NUMBERING = set()
    GM = t0.numbering.gathering
    for i in GM:
        LOCAL_NUMBERING.update(GM[i].full_vector)
    assert num_of_dofs_in_this_core == len(LOCAL_NUMBERING), "must be the case!"
    return 1


# noinspection PyUnboundLocalVariable
def test_Naive_Numbering_NO6_1trace():
    """"""
    if RANK == MASTER_RANK:
        load = random.randint(10, 199)
        print(f"--- [test_Naive_Numbering_NO6_1trace] @ FC-load = {load} ...... ", flush=True)
    else:
        load = None
    load = COMM.bcast(load, root=MASTER_RANK)
    FC = random_FormCaller_of_total_load_around(load, EDM_pool=('chaotic',))

    t1 = FC('1-t', numbering_parameters={'scheme_name': 'Naive'})

    GM_TEW = t1.numbering.trace_element_wise
    for i in GM_TEW:
        assert i in t1.mesh.trace.elements, "must be the case"
        assert len(set(GM_TEW.keys())) == t1.mesh.trace.elements.num, "must be the case"

    GM_TEW = COMM.gather(GM_TEW, root=MASTER_RANK)
    if RANK == MASTER_RANK:
        END = 0
        for i in range(t1.mesh.trace.elements.global_num):  # go through all (global) trace elements
            for gm_core in GM_TEW:
                if i in gm_core:
                    fv = gm_core[i].full_vector
                    assert fv[0] == END, f"must be the case."
                    NEW_END = fv[-1] + 1

            END = NEW_END


    num_of_dofs_in_this_core = t1.numbering.num_local_dofs
    LOCAL_NUMBERING = set()
    GM = t1.numbering.gathering
    for i in GM:
        LOCAL_NUMBERING.update(GM[i].full_vector)
    assert num_of_dofs_in_this_core == len(LOCAL_NUMBERING), "must be the case!"
    return 1


def test_Naive_Numbering_NO4_2trace():
    """"""
    if RANK == MASTER_RANK:
        print("--- [test_Naive_Numbering_NO4_2trace] ...... ", flush=True)

    mesh = MeshGenerator('crazy_periodic')([2, 2, 2], EDM='debug')
    space = SpaceInvoker('polynomials')([('Lobatto', 2), ('Lobatto', 3), ('Lobatto', 1)])
    FC = FormCaller(mesh, space)
    t2 = FC('2-t', numbering_parameters='Naive')

    GM = t2.numbering.gathering
    num_of_dofs_in_this_core = t2.numbering.num_local_dofs
    LOCAL_NUMBERING = set()
    benchmark = np.array(
        [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21],
         [3, 4, 5, 0, 1, 2, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37],
         [38, 39, 40, 41, 42, 43, 8, 9, 6, 7, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55],
         [41, 42, 43, 38, 39, 40, 24, 25, 22, 23, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67],
         [68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 16, 17, 18, 19, 20, 21, 10, 11, 12, 13, 14, 15],
         [71, 72, 73, 68, 69, 70, 78, 79, 80, 81, 32, 33, 34, 35, 36, 37, 26, 27, 28, 29, 30, 31],
         [82, 83, 84, 85, 86, 87, 76, 77, 74, 75, 50, 51, 52, 53, 54, 55, 44, 45, 46, 47, 48, 49],
         [85, 86, 87, 82, 83, 84, 80, 81, 78, 79, 62, 63, 64, 65, 66, 67, 56, 57, 58, 59, 60, 61]]
    )
    for i in GM:
        assert np.all(GM[i].full_vector == benchmark[i, :])
        LOCAL_NUMBERING.update(GM[i].full_vector)

    assert num_of_dofs_in_this_core == len(LOCAL_NUMBERING), "must be the case!"

    return 1


if __name__ == '__main__':
    # mpiexec -n 6 python _3dCSCG\TESTS\unittest_Naive_numbering.py
    # test_Naive_Numbering_NO4_2trace()
    test_Naive_Numbering_NO5_0trace()
    test_Naive_Numbering_NO6_1trace()
