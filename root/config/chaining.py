
# -*- coding: utf-8 -*-


from mpi4py import MPI
cOmm = MPI.COMM_WORLD
sIze: int = cOmm.Get_size()
rAnk: int = cOmm.Get_rank()





def CHAINING(method, *args, **kwargs):
    """Let all cores do thing in a sequence; so next core will not start unless the previous core has
    sent him a message to start.

    :param method: method or function to be executed.
    :param args: Each core needs to have a copy.
    :param kwargs: Each core needs to have a copy.
    :return:
    """
    if sIze == 1:
        return method(*args, **kwargs)
    else:
        if rAnk == 0:
            R0 = method(*args, **kwargs)
            cOmm.send(1, dest=1, tag=0)
            return R0
        elif rAnk < sIze-1: # intermediate cores
            _ = cOmm.recv(source=rAnk-1, tag=rAnk-1)
            Ri = method(*args, **kwargs)
            cOmm.send(1, dest=rAnk+1, tag=rAnk)
            return Ri
        elif rAnk == sIze-1: # the last core
            _ = cOmm.recv(source=rAnk-1, tag=rAnk-1)
            return method(*args, **kwargs)
        else:
            raise Exception()