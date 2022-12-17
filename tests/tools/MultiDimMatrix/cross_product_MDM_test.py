# -*- coding: utf-8 -*-
"""Here we use the MDM made from the cross product of standard forms to test
the MDM class.

@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/8/10 20:56
"""
import random
import sys

if './' not in sys.path: sys.path.append('./')

import numpy as np

from components.miscellaneous.miprint import miprint

from tests.objects.CSCG._2d.randObj.field import random_vector as rv2
from tests.objects.CSCG._2d.randObj.field import random_scalar as rs2
from tests.objects.CSCG._3d.randObj.field import random_vector as rv3

from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller

from objects.CSCG._2d.master import SpaceInvoker as space2
from objects.CSCG._2d.master import FormCaller as form2
from objects.CSCG._2d import mesh as mesh2

def test_MDM_sf_CrossProduct():
    """"""
    miprint("XXX [test_MDM_sf_CrossProduct] ...... ", flush=True)

    mesh = mesh2('crazy', c=0.0)([11, 12])
    space = space2('polynomials')([2, 2])

    fc2 = form2(mesh, space)

    s2 = rs2(fc2.mesh)
    v2 = rv2(fc2.mesh)
    V2 = rv2(fc2.mesh)

    w = fc2('0-f-o', hybrid=False)
    u = fc2('1-f-o', hybrid=False)
    v = fc2('1-f-o', hybrid=False)

    w.CF = s2
    w.CF.current_time = 0
    w.discretize()

    u.CF = v2
    u.CF.current_time = 0
    u.discretize()

    v.CF = V2
    v.CF.current_time = 0
    v.discretize()

    MDM = w.special.cross_product_1f__ip_1f(u, v, output='MDM')

    c = random.random()
    _2MDM = c * MDM
    __2MDM = MDM * c

    for _ in MDM:
        np.testing.assert_array_almost_equal(_2MDM[_], __2MDM[_])
        np.testing.assert_array_almost_equal(c*MDM[_], _2MDM[_])

    s2_X_v2__ip__V2 = s2.do.cross_product(v2).do.inner_product(V2)
    s2_X_v2__ip__V2.current_time = 0
    norm = s2_X_v2__ip__V2.do.compute_Ln_norm()
    NORM = MDM.do.evaluate([(w, w), (u, u), (v, v)])
    np.testing.assert_almost_equal(norm, NORM, decimal=4)

    mesh = MeshGenerator('crazy', c=0.0)([11, 12, 10])
    space = SpaceInvoker('polynomials')([('Lobatto', 2), ('Lobatto', 2), ('Lobatto', 2)])
    fc3 = FormCaller(mesh, space)

    v3 = rv3(fc3.mesh)
    V3 = rv3(fc3.mesh)
    s3 = rv3(fc3.mesh)

    w = fc3('1-f', hybrid=False)
    u = fc3('2-f', hybrid=False)
    v = fc3('2-f', hybrid=False)

    w.CF = s3
    w.CF.current_time = 0
    w.discretize()

    u.CF = v3
    u.CF.current_time = 0
    u.discretize()

    v.CF = V3
    v.CF.current_time = 0
    v.discretize()

    MDM = w.special.cross_product_2f__ip_2f(u, v, output='MDM')

    c = random.random()
    _2MDM = c * MDM
    __2MDM = MDM * c

    for _ in MDM:
        np.testing.assert_array_almost_equal(_2MDM[_], __2MDM[_])
        np.testing.assert_array_almost_equal(c*MDM[_], _2MDM[_])

    s2_X_v2__ip__V2 = s3.do.cross_product(v3).do.inner_product(V3)
    s2_X_v2__ip__V2.current_time = 0
    norm = s2_X_v2__ip__V2.do.compute_Ln_norm()
    NORM = MDM.do.evaluate([(w, w), (u, u), (v, v)])

    np.testing.assert_almost_equal(norm, NORM, decimal=3)

    return 1


if __name__ == '__main__':
    # mpiexec -n 4 python tests/tools/MultiDimMatrix/cross_product_MDM_test.py
    test_MDM_sf_CrossProduct()
