import numpy as np
#Find the index of the nearest value in an array
def find_nearest(array,values):
    indicies=np.zeros((len(values)),dtype=int)
    for i in range(len(values)): indicies[i] = np.abs(array-values[i]).argmin()
    return indicies
