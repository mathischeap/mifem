

"""


$ mpiexec -n 12 python TESTS\GENERAL\mpi4py_tests.py

"""

import sys
if './' not in sys.path: sys.path.append('../others/')



from root.config.main import *




a = np.array([1,1,1], dtype=float)




if RANK == MASTER_RANK:
    A = np.empty((3,), dtype=float)
    print(A)
else:
    A = None


COMM.Reduce(a, A, op=MPI.SUM, root=MASTER_RANK)


if RANK == MASTER_RANK:
    print(A)



b = np.array([2,2,2], dtype=float)
B = np.empty((3,), dtype=float)


if RANK == MASTER_RANK:
    COMM.Send([b, MPI.FLOAT], dest=SECRETARY_RANK, tag=1)

if RANK == SECRETARY_RANK:
    COMM.Recv([B, MPI.FLOAT], source=MASTER_RANK, tag=1)