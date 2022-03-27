
import pickle



def chain(filename):
    if filename[-3:] != '.mi': filename += '.mi'
    with open(filename, 'rb') as inputs:
        obj_dict = pickle.load(inputs)
    inputs.close()
    return obj_dict