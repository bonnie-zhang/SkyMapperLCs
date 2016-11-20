'''
from Anais Moller: https://github.com/anaismoller/Python_Utilities

Utility to perform a KDE with varying width according to uncertainties in measurements
Input:
    xlist= array xaxis where the KDE will be performed
    data_array= yaxis
    sigma_array= error array
    bandwidth=fixed initial bandwidth
Returns: 
    KDE array
'''

import numpy as np

def solve_gaussian(val,data_array,sigma_array,bandwidth):
    return (1. / np.sqrt(bandwidth**2 + sigma_array**2)) * np.exp(- (val - data_array) * (val - data_array) / (2 * (bandwidth**2 + sigma_array**2)))


def solve_kde(xlist,data_array,sigma_array,bandwidth):
    '''
    KDE with varying width according to unvertainties
    '''
    kde_array = np.array([])
    n = len(data_array)
    for xx in xlist:
        single_kde = solve_gaussian(xx,data_array,sigma_array,bandwidth)
        if np.ndim(kde_array) == 3:
            kde_array = np.concatenate((kde_array,single_kde[np.newaxis,np.newaxis,:]),axis=0)
        else:
            kde_array = np.dstack(single_kde)
    return kde_array
