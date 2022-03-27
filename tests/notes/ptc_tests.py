"""
Here we do some tests for the ptc project at mathischeap.com

# mpiexec -n 10 python TESTS\GENERAL\ptc_tests.py

"""
import sys
if './' not in sys.path: sys.path.append('../others/')

import numpy as np

from _3dCSCG.master import MeshGenerator, SpaceInvoker, FormCaller


# mesh = MeshGenerator('crazy', c=0.2, bounds=((0, 1), (0, 1), (0, 1)))([1, 1, 1])
#
# if 0 in mesh.elements:
#     e = mesh.elements[0]
# else:
#     exit()
#
# xi = np.linspace(-0.9, 0.6, 3)
# et = np.linspace(-0.7, 0.9, 3)
# sg = np.linspace(-0.8, 0.7, 3)
# xi, et, sg = np.meshgrid(xi, et, sg, indexing='ij')
#
# ct = e.coordinate_transformation
#
# x, y, z = ct.mapping(xi, et, sg)
#
#
# detJ = ct.Jacobian(xi, et, sg)
# # print(detJ)
#
#
# iG = ct.inverse_metric_matrix(xi, et, sg)
# # print(iG[1][2])
#
#
# iJ = ct.inverse_Jacobian_matrix(xi, et, sg)
# # print(iJ[1][2])
# # print(iJ[2][2])
#
# J = ct.Jacobian_matrix(xi, et, sg)
# # print(J[0][0])
#
#
# space = SpaceInvoker('polynomials')(
#     [('Lobatto', 2), ('Lobatto', 2), ('Lobatto', 2)])
# FC = FormCaller(mesh, space)
#
# f0 = FC('0-f', is_hybrid=False)
# f1 = FC('1-f', is_hybrid=False)
# f2 = FC('2-f', is_hybrid=False)
# f3 = FC('3-f', is_hybrid=False)
#
#
# M0 = f0.operators.inner(f0, quad_degree=(4,4,4))
# print(M0[0][13,:].toarray())
#
#
# M1 = f1.operators.inner(f1, quad_degree=(4,4,4))
# print(M1[0][25,:].toarray())
#
#
#
# M2 = f2.operators.inner(f2, quad_degree=(4,4,4))
# print(M2[0][30,:].toarray())
#
#
#
# M3 = f3.operators.inner(f3, quad_degree=(4,4,4))
# print(M3[0].toarray())






# generate a random unit vector

n = np.random.rand(3)
n = n / np.sqrt(np.sum(n**2))


u = np.random.rand(3)





u_dot_n = np.dot(u, n)

perp = u_dot_n * n

u_cross_n = np.cross(u, n)


tang = np.cross(n, u_cross_n)


print(tang + perp - u)


print(np.dot(tang, n))


print(np.dot(u_cross_n, n))





