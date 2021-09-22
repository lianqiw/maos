#!/usr/bin/env python3
import glob
import os
import numpy as np
np.set_printoptions(threshold=100,suppress=False,precision=4,floatmode='maxprec',linewidth=120)
from scipy.special import erf
from numpy import sqrt, exp, log, floor, ceil, nan, pi
from numpy.random import rand, randn
import matplotlib.pyplot as plt
try:
    from natsort import natsorted
except:
    natsorted=sorted
try:
    from libaos import *
except:
    import lib2py #it generates libaos
    from libaos import *

from readbin import readbin
from draw import draw, locembed
import maos_client

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

def cummean(y,axis=0):
    ys=y.shape
    x=1/np.arange(1,1+ys[axis])
    if len(ys)==2:
        if axis==0:
            x.shape=(x.shape[0],1)
        else:
            x.shape=(1,x.shape[0])
    
    yc=np.cumsum(y,axis=axis)
    return yc*x
def maos_cumu(files, seeds=None, nsim0=0): ##return cumulative average
    res,fds=maos_res(files,seeds,0,0)
    print(fds)
    nsim=res.shape[-1]
    if nsim0<=0:
        nsim0=min(5000, np.int(nsim*0.1))
    yy=np.arange(1, nsim+1-nsim0)
    xx=nsim0+yy
    yy.shape=(1,1,nsim-nsim0)
    resc=np.cumsum(res[:,:,nsim0:], axis=2)/yy
    resc[res[:,:,nsim0:]==0]=nan
    return resc,xx
def maos_res(fds, seeds=None, iframe1=0.2, iframe2=1):
    return maos_res_do(fds, "Res", seeds, iframe1, iframe2)
def maos_res_each_old(fds, seeds=None, iframe1=0.2, iframe2=1):
    return maos_res_do(fds, "Resclep", seeds, iframe1, iframe2)
def maos_res_each(fds, seeds=None, iframe1=0.2, iframe2=1):
    return maos_res_do(fds, "extra/Resp", seeds, iframe1, iframe2)
def maos_res_do(fdin, name, seeds=None, iframe1=0.2, iframe2=1):
    if type(fdin) is list:
        fds2=[]
        for fdini in fdin:
            fds2+=natsorted(glob.glob(fdini+"/",recursive=1))
    else:
        fds2=natsorted(glob.glob(fdin+"/",recursive=1))
    
    fds=[]
    resall=None
    for fd in fds2: #loop over directory
        if seeds is None:
            fns=natsorted(glob.glob(fd+"/"+name+"_*.bin"))
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
            res=read(fn)
            if res is None or res.shape[0]==0:
                continue
            if name=="Res":
                if res[3].size>0: #split tomography
                    res=res[3]
                    split=1
                else: #integrated
                    res=res[2]
                    split=0
            elif name=="Resclep":
                res=np.stack(res, axis=1)
                split=0
            elif name=="extra/Resp":
                res=np.stack(res[-1],axis=1)
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
            if n1 < n2 and n2 <= res.shape[0]:
                res=res[n1:n2]
                res=res.mean(0,keepdims=1)
            else:
                res=res.T
                res.shape=(1,res.shape[0],res.shape[1])
            #if max(res[-1]).all()==0:
            #    print(fn, 'Incomplete result')
            #    continue #unfinished result
    
            mres+=res*1e18
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
    print(len(fds), ' results are read')
    return (resall,np.array(fds))
    #else:
    #return resall
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

"""
import inspect, dis
def nargout():
    "Return how many values the caller is expecting"
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
"""
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

#remove tip/tilt/focus from gradients
def grad_ttfr(saloc,g):
    if g.dtype==object:
        gv=np.empty(g.shape, dtype=object)
        print(gv.shape)
        for ig,gi in np.ndenumerate(g):
            gv[ig]=grad_ttfr(saloc,gi)
        return gv
    if saloc.shape[0]==2: #2*nsa
        nsa=saloc.shape[1]
        tt=saloc.flatten() 
    elif saloc.shape[1]==2: #nsa*2
        nsa=saloc.shape[0]
        tt=saloc.T.flatten()
    else:
        raise(ValueError('saloc should be 2xnsa or nsax2'))
    
    if g.shape[0]==2: #2*nsa
        gv=g.flatten('C')
    elif g.shape[0]==2*nsa:
        gv=g
    elif g.shape[1]==2:
        gv=g.flatten('F')
    else:
        raise(ValueError('g should bd 2*nsa, nsa*2 or nsa*2*m'))
    
    mod=np.c_[np.ones((nsa*2,1)), tt] #nsa*3
    
    rmod=np.linalg.pinv(mod)
    
    ptt=rmod@gv
    g2v=gv-mod@ptt
    g2v.shape=[2,nsa]
    return g2v
#remove zernike modes from rmin to rmax from 1-D OPD and loc
def opd_loc_project_zernike(opd, loc, rmin, rmax, radonly=0):
    D=np.max(np.max(loc,axis=1)-np.min(loc, axis=1))
    mod=zernike(loc, D, rmin, rmax, radonly).T
    rmod=np.linalg.pinv(mod)
    return mod@(rmod@opd)
    
def opd_loc_remove_zernike(opd, loc, rmin, rmax, radonly=0):
    D=np.max(np.max(loc,axis=1)-np.min(loc, axis=1))
    mod=zernike(loc, D, rmin, rmax, radonly).T
    rmod=np.linalg.pinv(mod)
    opd2=opd-mod@(rmod@opd)
    return opd2
#remove zernike modes from rmin to rmax from 2-D OPD
def opd_remove_zernike(opd, mask, rmin, rmax, radonly=0):
    '''OPD is 2d array, mask indicates valid points, rmin and rmax are minimum and maximum zernike order (inclusive)'''

    if mask is None:
        mask=opd!=0
    dx=0.4 #doesn't matter
    oloc=mksqloc(opd.shape[1], opd.shape[0], dx, dx, -opd.shape[1]/2*dx, -opd.shape[0]/2*dx)
    oloc.shape=[2,opd.shape[0],opd.shape[1]]
    oloc2=oloc[:,mask].copy() #copy is necessary for C to access the data
    opd2=np.zeros(opd.shape)
    opd2[mask]=opd_loc_remove_zernike(opd[mask], oloc2, rmin, rmax, radonly);
    #D=np.max(np.max(oloc2,axis=1)-np.min(oloc2, axis=1))
    #mod=zernike(oloc2, D, rmin, rmax, radonly).T
    #rmod=np.linalg.pinv(mod)
    #opd2[mask]=opd[mask]-mod@(rmod@opd[mask])
    return opd2
def opd_loc_remove_focus(opd, loc):
    return opd_loc_remove_zernike(2,2,1)
def read_many(fdin):
    fds2=natsorted(glob.glob(fdin,recursive=1))
    fds=[]
    res=[]
    for fd in fds2: 
        try:
            tmp=read(fd)
            fds.append(fd)
            res.append(tmp)
        except:
            pass
    return np.array(res),fds
def eye(nx, val=1):
    dd=np.zeros((nx,nx))
    np.fill_diagonal(dd, val)
    return dd
def svd_inv(A, thres=0, tikcr=0):
    '''svd_inv(A, thres=0, tikcr=0)'''
    u,s,v=np.linalg.svd(A)
    if tikcr>0:
        dd=np.zeros(A.shape)
        np.fill_diagonal(dd,s[0]*tikcr)
        u,s,v=np.linalg.svd(A+eye(A.shape[0], s[0]*tikcr))
    si=1./s
    if thres>0:
        si[s<s[0]*thres]=0
    return (u*si)@v