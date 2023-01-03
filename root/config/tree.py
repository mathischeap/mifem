# -*- coding: utf-8 -*-
from mpi4py import MPI
COMM = MPI.COMM_WORLD
SIZE: int = COMM.Get_size()
RANK: int = COMM.Get_rank()
import numpy as np


def TREE(factor=2):
    """We distribute our cores with a tree structure of branch ``factor``.

    :param int factor: (`default`: ``2``).
    :return: A generator.
    """
    if factor == 2:
        H = np.array([(i, i+1) if i + 1 < SIZE else (i, -1) for i in range(0, SIZE, 2)])
        yield ___parse_senders_and_receivers_factor2___(H)
        KEYS = H[:, 0]
        LEN = np.size(KEYS)
        while LEN > 1:
            H = np.array([(KEYS[i], KEYS[i + 1]) if i + 1 < LEN else (KEYS[i], -1) for i in range(0, LEN, 2)])
            yield ___parse_senders_and_receivers_factor2___(H)
            KEYS = H[:, 0]
            LEN = np.size(KEYS)
    else:
        raise Exception(f"factor={factor} not coded.")


def ___parse_senders_and_receivers_factor2___(H):
    """"""
    for m, row in enumerate(H):
        if row[1] == -1:
            pass
        else:
            i, j = row
            i = int(i)
            j = int(j)
            tag = m
            if RANK == i:
                recv_send = ('recv', {'source': j, 'tag': tag})
                return recv_send
            if RANK == j:
                recv_send = ('send', {'dest': i, 'tag': tag})
                return recv_send
