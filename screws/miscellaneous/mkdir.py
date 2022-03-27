"""make a dir"""


from root.config.main import rAnk, mAster_rank
import os

def mkdir(folder_name):
    """"""
    if rAnk == mAster_rank:
        if os.path.isdir(folder_name):
            pass
        else:
            os.mkdir(folder_name)