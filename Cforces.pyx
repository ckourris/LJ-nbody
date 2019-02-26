# This library is a python wrapper to calculate get_forces using
# C(++) accelerated code.

cimport numpy as np
np.import_array()

# Include python cfunction
cdef extern from "getforces.h":
    void getforces(double* inarray, double* out_array, int N, double boxdim, double cutoff)

# Python wrapper for cfunction
def c_getforces(np.ndarray[double, ndim=2, mode='c'] in_array not None,
                np.ndarray[double, ndim=2, mode='c'] out_array not None,
                boxdim, cutoff):
    getforces(<double*> np.PyArray_DATA(in_array), <double*> np.PyArray_DATA(out_array),
              in_array.shape[0], boxdim, cutoff);
    return
