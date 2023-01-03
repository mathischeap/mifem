# -*- coding: utf-8 -*-
from mpi4py import MPI
COMM = MPI.COMM_WORLD
SIZE: int = COMM.Get_size()
RANK: int = COMM.Get_rank()


def CHAINING(method, *args, **kwargs):
    """Let all cores do thing in a sequence; so next core will not start unless the previous core has
    sent him a message to start.

    :param method: method or function to be executed.
    :param args: Each core needs to have a copy.
    :param kwargs: Each core needs to have a copy.
    :return:
    """
    if SIZE == 1:
        return method(*args, **kwargs)
    else:
        if RANK == 0:
            R0 = method(*args, **kwargs)
            COMM.send(1, dest=1, tag=0)
            return R0
        elif RANK < SIZE-1:  # intermediate cores
            _ = COMM.recv(source=RANK - 1, tag=RANK - 1)
            Ri = method(*args, **kwargs)
            COMM.send(1, dest=RANK + 1, tag=RANK)
            return Ri
        elif RANK == SIZE-1:  # the last core
            _ = COMM.recv(source=RANK - 1, tag=RANK - 1)
            return method(*args, **kwargs)
        else:
            raise Exception()
