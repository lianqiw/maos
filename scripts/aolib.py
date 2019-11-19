#!/usr/bin/env python3
try:
    from libaos import *
except:
    import lib2py #it generates libaos
    from libaos import *
import glob
import os
import numpy as np
np.set_printoptions(threshold=100,suppress=False,precision=4,floatmode='maxprec',linewidth=120)


import matplotlib as mpl
from cycler import cycler
mpl.rcParams['axes.prop_cycle'] = cycler(color='bgrcmyk')
mpl.rcParams['axes.grid']=True
#mpl.rcParams['grid.color'] = 'k'
#mpl.rcParams['grid.linestyle'] = '--'
#mpl.rcParams['grid.linewidth'] = 0.5

mpl.rcParams['lines.linewidth'] = 1.0
mpl.rcParams['lines.dashed_pattern'] = [6, 4]
mpl.rcParams['lines.dashdot_pattern'] = [6, 3, 1, 3] #dash off dot off
mpl.rcParams['lines.dotted_pattern'] = [1, 3]
mpl.rcParams['lines.scale_dashes'] = False

mpl.rcParams['axes.xmargin']=0
mpl.rcParams['axes.ymargin']=0
mpl.rcParams['axes.autolimit_mode']='round_numbers'

mpl.rcParams['font.size']=10
mpl.rcParams['savefig.dpi']=120
mpl.rcParams['image.cmap']='jet'

mpl.rcParams['figure.autolayout']=True
#mpl.rcParams['figure.subplot.left']=0.1
#mpl.rcParams['figure.subplot.right']=0.9

#For ploting
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, semilogx, semilogy, loglog, xlabel, ylabel, legend, grid, clf, figure, subplot, xlabel, ylabel, title, xlim, ylim, close, savefig
from scipy.special import erf
from numpy import sqrt, exp, log, floor, ceil, nan
from numpy.random import rand, randn

plt.ion() #enable interactive mode.


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
    fds2=glob.glob(fds)
    resall=None
    for fd in fds2:
        if seeds is None:
            fns=glob.glob(fd+"/Res_*.bin")
        else:
            fns=list()
            if type(seeds)!=list:
                seeds=[seeds]
            for seed in seeds:
                fns.append(fd+'/Res_{}.bin'.format(seed))
        nseed=0
        mres=0
        split=-1
        for fn in fns:
            if not os.path.exists(fn):
                print(fn, 'does not exist')
                continue
            res=readbin(fn)
            if res is None:
                continue
            if res[3].size>0: #split tomography
                res=res[3]
                split=1
            else: #integrated
                res=res[2]
                split=0
            
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
            mres=mres*(1/nseed)
            if resall is None:
                resall=mres
            else:
                resall=np.vstack((resall, mres))
    if resall is None:
        resall=np.array([nan,nan,nan])
    if len(fds2)>1:
        return (resall,fds2)
    else:
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
    lines=['-','--','-.',':'];
    return lines[np.mod(ii,len(lines))]

def reset_color():
    plt.gca().set_prop_cycle(None)
