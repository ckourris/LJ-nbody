# This library is a python wrapper to calculate get_forces using
# C(++) accelerated code.

cimport numpy as np
np.import_array()

# Include python cfunction
cdef extern from "accelerate.h":
    void getforces(double* inarray, double* out_array, int N, double boxdim, double cutoff)
    void getenergies(double* inarray_v, double* inarray_pos, double* out_array,
                    int N, double boxdim, double cutoff)

# Python wrapper for cfunction
def c_getforces(np.ndarray[double, ndim=2, mode='c'] in_array not None,
                np.ndarray[double, ndim=2, mode='c'] out_array not None,
                boxdim, cutoff):
    getforces(<double*> np.PyArray_DATA(in_array), <double*> np.PyArray_DATA(out_array),
              in_array.shape[0], boxdim, cutoff);
    return

def c_getenergies(np.ndarray[double, ndim=2, mode='c'] v_array not None,
                  np.ndarray[double, ndim=2, mode='c'] pos_array not None,
                  np.ndarray[double, ndim=1, mode='c'] out_array not None,
                  boxdim, cutoff):
    getenergies(<double*> np.PyArray_DATA(v_array), <double*> np.PyArray_DATA(pos_array),
              <double*> np.PyArray_DATA(out_array), v_array.shape[0], boxdim, cutoff);
    return
