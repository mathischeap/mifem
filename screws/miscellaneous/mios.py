# -*- coding: utf-8 -*-
""""""


from root.config.main import rAnk, mAster_rank
import os

def mkdir(folder_name):
    """"""
    if rAnk == mAster_rank:
        if os.path.isdir(folder_name):
            pass
        else:
            os.mkdir(folder_name)

def remove(*file_names):
    if rAnk == mAster_rank:
        for file_name in file_names:
            os.remove(file_name)


def rmdir(folder_name):
    if rAnk == mAster_rank:
        os.rmdir(folder_name)