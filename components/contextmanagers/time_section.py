# -*- coding: utf-8 -*-
""""""

import contextlib
from time import localtime, strftime, time
from root.config.main import *


@contextlib.contextmanager
def timesection(info=None):
    COMM.barrier()
    if RANK == MASTER_RANK:
        print(" <{}> starts at [".format(info) + strftime("%Y-%m-%d %H:%M:%S", localtime()) + ']')
    ts = time()
    yield
    COMM.barrier()
    if RANK == MASTER_RANK:
        print(" <{}> ends at [".format(info) + strftime("%Y-%m-%d %H:%M:%S", localtime()) + ']')
        print(" <{}> costs: [%.5f seconds]\n".format(info) % (time()-ts), flush=True)
