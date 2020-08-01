#!/usr/bin/env python3
import glob
import os
import sys
import numpy as np
np.set_printoptions(threshold=100,suppress=False,precision=4,floatmode='maxprec',linewidth=120)
from scipy.special import erf
from numpy import sqrt, exp, log, floor, ceil, nan
from numpy.random import rand, randn

try:
    from libaos import *
except:
    import lib2py #it generates libaos
    from libaos import *

from readbin import readbin

from draw import draw


#To dock multiple figures. Does not work very well.
def dock_figure():
    #mpl.use('module://mpldock.backend')
    from mpldock import persist_layout
    persist_layout('maos')
    plt.switch_backend('module://mpldock.backend')
    ipython.magic('gui qt5') #replacing plt.show() that blocks


from matplotlib.pyplot import plot, semilogx, semilogy, loglog, xlabel, ylabel, legend, grid, clf, subplot, xlabel, ylabel, title, xlim, ylim, close, savefig

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

def maos_res(fds, seeds=None, iframe1=0.2, iframe2=1):
    return maos_res_do(fds, "Res", seeds, iframe1, iframe2)
def maos_res_each(fds, seeds=None, iframe1=0.2, iframe2=1):
    return maos_res_do(fds, "Resclep", seeds, iframe1, iframe2)
def maos_res_do(fdin, name, seeds=None, iframe1=0.2, iframe2=1):
    fds2=glob.glob(fdin+"/",recursive=1)
    fds=[]
    resall=None
    for fd in fds2: #loop over directory
        if seeds is None:
            fns=glob.glob(fd+"/"+name+"_*.bin")
        else:
            fns=list()
            for seed in seeds:
                fns.append(fd+'/{}_{}.bin'.format(name,seed))

        nseed=0
        mres=0
        split=-1
        for fn in fns: #loop over seed
            if not os.path.exists(fn):
                print(fn, 'does not exist')
                continue
            res=readbin(fn)
            if res is None:
                continue
            if len(res)>4: #per direction, take on axis
                res=res[0]
                split=0
            elif len(res)==4:
                if res[3].size>0: #split tomography
                    res=res[3]
                    split=1
                else: #integrated
                    res=res[2]
                    split=0
            else:
                print('Invalid result')
                continue
            if iframe1<1:
                n1=round(iframe1*res.shape[0])
            else:
                n1=iframe1
            if iframe2<=1:
                n2=round(iframe2*res.shape[0])
            else:
                n2=iframe2
            res=res[n1:n2,:]
            if max(res[-1,:])==0:
                print(fn, 'Incomplete result')
                continue #unfinished result
            res=res.mean(0)*1e18
            mres+=res
            nseed+=1
        if nseed>0:
            fds.append(fd[0:-1])
            mres=mres*(1/nseed)
            if resall is None:
                resall=mres
            else:
                resall=np.vstack((resall, mres))
        else:
            print(fd, fns, nseed, ' has no valid results')
    if resall is None:
        resall=np.array([nan,nan,nan])
        if not os.path.exists(fdin):
            print(fdin, 'does not exist')
        else:
            print(fdin, ' has no valid results')
    #if len(fds)>1:
    #    print(*fds, sep="\n")
    return (resall,fds)
    #else:
    return resall
def mysqrt(x):
    if type(x) is np.ndarray:
        return np.sign(x)*np.sqrt(np.abs(x))
    elif x<0:
        return -np.sqrt(-x)
    else:
        return np.sqrt(x)

def rss(*args):
    arr=np.array(args)
    return mysqrt(np.sum(np.sign(arr)* arr**2))
    
def styles(ii):
    reset_color()
    lines=['-','--','-.',':'];
    return lines[np.mod(ii,len(lines))]

def reset_color():
    plt.gca().set_prop_cycle(None)


import inspect, dis
def nargout():
    """Return how many values the caller is expecting"""
    f = inspect.currentframe()
    f = f.f_back.f_back
    c = f.f_code
    i = f.f_lasti
    bytecode = c.co_code
    instruction = bytecode[i+3]
    if instruction == dis.opmap['UNPACK_SEQUENCE']:
        howmany = bytecode[i+4]
        return howmany
    elif instruction == dis.opmap['POP_TOP']:
        return 0
    return 1

def figure(*args, **kargs):
    plt.figure(*args, **kargs)
    plt.clf()

    
def width_at(x, y, height):
    return len(np.argwhere(y>=max(y)*height))*(x[1]-x[0])

def cellsum(x):
    xsum=np.zeros(x.shape)
    for ii in range(x.size):
        jj=np.unravel_index(ii, x.shape)
        xsum[jj]=x[jj].sum()

    return xsum
