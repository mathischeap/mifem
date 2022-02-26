

"""


$ mpiexec -n 12 python TESTS\GENERAL\mpi4py_tests.py

"""

import sys
if './' not in sys.path: sys.path.append('./')



from root.config import *




a = np.array([1,1,1], dtype=float)




if rAnk == mAster_rank:
    A = np.empty((3,), dtype=float)
    print(A)
else:
    A = None


cOmm.Reduce(a, A, op=MPI.SUM, root=mAster_rank)


if rAnk == mAster_rank:
    print(A)



b = np.array([2,2,2], dtype=float)
B = np.empty((3,), dtype=float)


if rAnk == mAster_rank:
    cOmm.Send([b, MPI.FLOAT], dest=sEcretary_rank, tag=1)

if rAnk == sEcretary_rank:
    cOmm.Recv([B, MPI.FLOAT], source=mAster_rank, tag=1)