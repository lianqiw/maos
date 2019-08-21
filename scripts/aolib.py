#!/usr/bin/env python3
try:
    from libaos import *
except:
    import lib2py #it generates libaos
    from libaos import *

def iscell(arr):
    if type(arr)==np.ndarray and arr.dtype.name=='object':
        return True
    else:
        return False

'''use of .flat with an FORTRAN ordered array will lead to non-optimal memory
access as adjacent elements in the flattened array (iterator, actually) are not
contiguous in memory'''

def isequal(a, b):
    if type(a)==np.ndarray and type(b)==np.ndarray:
        if a.dtype.name=='object':
            for ii in range(a.size):
                if not isequal(a.item(ii),b.item(ii)):
                    return False
            return True
        else:
            an=np.linalg.norm(a)
            bn=np.linalg.norm(b)
            return abs(an-bn)<(an+bn)*1e-15
    elif sp.issparse(a) and sp.issparse(b):
        x=np.random.rand(a.shape[1])
        return isequal(a*x, b*x)
    else:
        return a==b


def test_read():
    import matplotlib.pyplot as plt
    import os
    path=b'/home/lianqiw/work/aos/comp/optim/bin/test/setup/'
    
    loc=readbin(path+b'aloc.bin')[0,0]
    opd=np.arange(0,loc.shape[1])
    opdmap=loc_embed(loc, opd);
    plt.imshow(opdmap)

    #path=b'/home/lianqiw/.aos/cache/SELOTF/'
    
    
    for file in os.listdir(path):
        if file.endswith(b"_py.bin") or not file.endswith(b".bin"):
            continue


        tmp=readbin(path+file)
        file2=file[0:-4]+b'_py.bin'
        if not os.path.isfile(path+file2):
            writebin(tmp, path+file2)
        tmp2=readbin(path+file2)
        if isequal(tmp, tmp2):
            print(file, 'equal')
        else:
            raise(Exception(file, 'not equal'))
    
def test_mkdtf():
    out=mkdtf([0.5e-6], 1./64., 2., 64, 64, 8,8,10e-6,10e-6,[],[],0,[],0,0)
    return out

if __name__ == '__main__':
    a=readbin('non.bin');
    b=test_mkdtf()
    print(a)
    print(b)
