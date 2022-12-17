# -*- coding: utf-8 -*-
"""Here, we store some most fundamental functions for mifem.

In this script, we DO NOT use the structure of naming files and folders of the mifem library.

"""
from components.miscellaneous.timer import check_filename_mi
from root.config.main import MASTER_RANK, RANK, COMM
import pickle
from root.read.main import read


def save(obj, filename):
    """

    Parameters
    ----------
    obj
    filename

    Returns
    -------

    """
    filename = check_filename_mi(filename)
    if isinstance(obj, (list, tuple)):
        _2bs_ = list()
        _sif_ = tuple()
        for obj_i in obj:
            bs = obj_i.___PRIVATE_save___(filename, do_save=False)
            _2bs_.append(bs)

            sif = obj_i.___PRIVATE_saving_info___()
            _sif_ += (sif,)

        _2bs_.append(_sif_)

        COMM.barrier()
        if RANK == MASTER_RANK:
            with open(filename, 'wb') as output:
                pickle.dump(_2bs_, output, pickle.HIGHEST_PROTOCOL)
            output.close()
    else:
        obj.___PRIVATE_save___(filename, do_save=True)



if __name__ == '__main__':
    rd = read
