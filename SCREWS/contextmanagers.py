# -*- coding: utf-8 -*-



import contextlib
from time import localtime, strftime, time
from root.config import *



@contextlib.contextmanager
def timesection(info=None):
    cOmm.barrier()
    if rAnk == mAster_rank:
        print(" <{}> starts at [".format(info) + strftime("%Y-%m-%d %H:%M:%S", localtime()) + ']')
    ts = time()
    yield
    cOmm.barrier()
    if rAnk == mAster_rank:
        print(" <{}> ends at [".format(info) + strftime("%Y-%m-%d %H:%M:%S", localtime()) + ']')
        print(" <{}> costs: [%.5f seconds]\n".format(info)%(time()-ts), flush=True)