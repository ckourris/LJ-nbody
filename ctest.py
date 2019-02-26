from cforce import *
import numpy as np
from Utilities import *
from Particle3D import*

def pforce(a,b,boxdim, cutoff):
    at = Particle3D(pos=a)
    bt = Particle3D(pos=b)
    return LJ_Force(Particle3D.pbc_sep(at, bt, boxdim), cutoff)

N=1000

testdata = np.random.rand(N,2,3)*10-5
testdim = np.random.random(N)*5
testcutoff = np.random.random(N)*10

for i in range(N):
    a = testdata[i,0]
    b = testdata[i,1]
    #print('Test: ', 1)
    #print('a = ', a, 'b = ', b)
    pf = pforce(a,b,testdim[i],testcutoff[i])
    #print('Python force: ', pf)
    #print('C++: ')
    c = np.zeros((2,3), dtype=np.float64)
    c_getforces(testdata[i], c, testdim[i], testcutoff[i])
    #print(c)
    if(np.linalg.norm(pf-c[0]) > 0.0001):
        print("ERRRRRRRRROOOOOORRRRRRRRR", i)
        print(a, b, testdim[i], testcutoff[i])
        print(pf)
        print(c)
