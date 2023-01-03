# -*- coding: utf-8 -*-
""""""
from root.config.main import RANK, MASTER_RANK, COMM
import os


def isfile(filename):
    """"""
    if RANK == MASTER_RANK:
        ToF = os.path.isfile(filename)
    else:
        ToF = None

    return COMM.bcast(ToF, root=MASTER_RANK)


def mkdir(folder_name):
    """"""
    if RANK == MASTER_RANK:
        if os.path.isdir(folder_name):
            pass
        else:
            os.mkdir(folder_name)


def remove(*file_names):
    COMM.barrier()
    if RANK == MASTER_RANK:
        for file_name in file_names:
            os.remove(file_name)


def rmdir(folder_name):
    if RANK == MASTER_RANK:
        os.rmdir(folder_name)


def listdir(folder_name):
    """Return all filenames in the folder."""
    if RANK == MASTER_RANK:
        return os.listdir(folder_name)


def cleandir(folder_name):
    """clean all files in a folder.

    Use this function very carefully. You do not want to delete all your important files accidentally.

    """
    COMM.barrier()
    if RANK == MASTER_RANK:
        files = listdir(folder_name)
        for file in files:
            os.remove(folder_name + '/' + file)
