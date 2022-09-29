# -*- coding: utf-8 -*-
""""""


from root.config.main import rAnk, mAster_rank, cOmm
import os

def mkdir(folder_name):
    """"""
    if rAnk == mAster_rank:
        if os.path.isdir(folder_name):
            pass
        else:
            os.mkdir(folder_name)

def remove(*file_names):
    cOmm.barrier()
    if rAnk == mAster_rank:
        for file_name in file_names:
            os.remove(file_name)


def rmdir(folder_name):
    if rAnk == mAster_rank:
        os.rmdir(folder_name)

def listdir(folder_name):
    """Return all filenames in the folder."""
    if rAnk == mAster_rank:
        return os.listdir(folder_name)

def cleandir(folder_name):
    """clean all files in a folder.

    Use this function very carefully. You do not want to delete all your important files accidentally.

    """
    cOmm.barrier()
    if rAnk == mAster_rank:
        files = listdir(folder_name)
        for file in files:
            os.remove(folder_name + '/' + file)