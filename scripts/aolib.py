#!/usr/bin/env python3
import glob
import os
import numpy as np
np.set_printoptions(threshold=100,suppress=False,precision=4,floatmode='maxprec',linewidth=120)
from scipy.special import erf
from scipy.integrate import cumtrapz
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

from readbin import readbin, headers
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

def cummean(y,axis=0,skip=0,rms=0):
    ys=y.shape
    x=1/np.arange(1,1+ys[axis]-skip)
    yc=np.zeros(y.shape)
    if len(ys)==2:
        if axis==0:
            x.shape=(x.shape[0],1)
            y2=y[skip:]
            yc2=yc[skip:]
        else:
            x.shape=(1,x.shape[0])
            y2=y[:,skip:]
            yc2=yc[skip:,:]
    else:
        y2=y[skip:]
        yc2=yc[skip:]
    if rms!=0:
        yc2[:]=np.sqrt(np.cumsum(y2**2,axis=axis)*x)
    else:
        yc2[:]=np.cumsum(y2,axis=axis)*x
    return yc
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
    '''Results are in order of High, T/T, NGS total, Focus'''
    return maos_res_do(fds, "Res", seeds, iframe1, iframe2)
def maos_res_tot(fds, seeds=None, iframe1=0.2, iframe2=1):
    '''Results are High + NGS tot'''
    res,fds=maos_res_do(fds, "Res", seeds, iframe1, iframe2)
    res2=res[:,0]+res[:,2]
    return res2,fds
def maos_res_hi(fds, seeds=None, iframe1=0.2, iframe2=1):
    '''Results are High order only'''
    res,fds=maos_res_do(fds, "Res", seeds, iframe1, iframe2)
    res2=res[:,0]
    return res2,fds
def maos_res_each_old(fds, seeds=None, iframe1=0.2, iframe2=1):
    return maos_res_do(fds, "Resclep", seeds, iframe1, iframe2)
def maos_res_each(fds, seeds=None, iframe1=0.2, iframe2=1):
    return maos_res_do(fds, "extra/Resp", seeds, iframe1, iframe2)
def maos_res_do(fdin, name, seeds=None, iframe1=0.2, iframe2=1):
    if type(fdin) is not list:
        fdin=[fdin]
    fds2=[]
    for fdini in fdin:
        fds2+=glob.glob(fdini+'/',recursive=1) #find only directory

    #natsorted not work well if there is trailing /
    fds2=natsorted([fd[0:-1] for fd in fds2])
    fds=[]
    resall=None
    nseed=0
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
            if type(res) is tuple:
                res=res[0]
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
            n_valid=np.nonzero(res[:,0]>0)[0][-1]                
            if iframe1<1:
                n1=round(iframe1*n_valid)
            else:
                n1=iframe1

            if iframe2<=1:
                n2=round(iframe2*n_valid)
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
            fds.append(fd)
            mres=mres*(1/nseed)
            if resall is None:
                resall=mres
            else:
                resall=np.vstack((resall, mres))
        else:
            print(fd, fns, nseed, ' has no valid results')
    if resall is None:
        resall=np.array([nan,nan,nan])
        #if not os.path.exists(fdin):
        #    print(fdin, 'does not exist')
        #else:
        print(fdin, ' has no valid results')
    #if len(fds)>1:
    #    print(*fds, sep="\n")
    print('{} results are read for {} seeds.'.format(len(fds), nseed))
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
def rms(*args):
    arr=np.array(args)
    return np.sqrt(np.mean(arr**2))
    
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
def grad_ttfr(grad, saloc):
    if grad.dtype==object:
        gv=np.empty(grad.shape, dtype=object)
        print(gv.shape)
        for ig,gi in np.ndenumerate(grad):
            gv[ig]=grad_ttfr(saloc,gi)
        return gv
    if saloc.shape[1]==2 and saloc.shape[0]!=2: #nsa*2
        saloc=saloc.T
        
    if saloc.shape[0]==2: #2*nsa
        nsa=saloc.shape[1]
        tt=saloc.flatten() 
    else:
        raise(ValueError('saloc should be 2xnsa or nsax2'))
    
    if grad.shape[0]==2: #2*nsa
        gv=grad.flatten('C')
    elif grad.shape[0]==2*nsa:
        gv=grad
    elif grad.shape[1]==2:
        gv=grad.flatten('F')
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
            tmp=readbin(fd)
            fds.append(fd)
            res.append(tmp)
        except:
            print('Fail to read',fd)
            pass
    return simplify(np.array(res)),fds
def read_many_dict(fdin):
    fds2=natsorted(glob.glob(fdin,recursive=1))
    res={}
    for fd in fds2: 
        try:
            res[fd]=readbin(fd)
        except:
            print('Fail to read',fd)
            pass
    return res
def eye(nx, val=1):
    dd=np.zeros((nx,nx))
    np.fill_diagonal(dd, val)
    return dd
def svd_inv(A, thres=0, tikcr=0):
    '''svd_inv(A, thres=0, tikcr=0)'''
    u,s,v=np.linalg.svd(A)
    if tikcr>0:
        u,s,v=np.linalg.svd(A+eye(A.shape[0], s[0]*tikcr))
    si=1./s
    if thres>0:
        si[s<s[0]*thres]=0
    return (u*si)@v
def svd_pinv(mod, thres=0, tikcr=0):
    '''svd_pinv(mod, thres=0, tikcr=0)'''
    '''Each row of mod is treated as a mod'''
    if mod.shape[1]<mod.shape[0]:
        print('mod has wrong shape, use its transpose instead:', mod.shape)
        return svd_pinv(mod.T, thres, tikcr)

    mmt=mod@mod.T 
    immt=svd_inv(mmt, thres, tikcr)
    return mod.T@immt

def plot_smooth(x,y,color=''):
    from scipy.interpolate import make_interp_spline, BSpline

    #define x as 200 equally spaced values between the min and max of original x 
    xnew = np.linspace(x.min(), x.max(), 200) 

    #define spline
    spl = make_interp_spline(x, y, k=3)
    y_smooth = spl(xnew)

    #create smooth line chart 
    plt.plot(xnew, y_smooth,color)

def radial_profile(data, center=None, enclosed=0):
    '''Compute the radial average or radially enclosed energy. radial start with 1'''
    if center is None:
        center=(data.shape[0]>>1, data.shape[1]>>1)
    y,x = np.indices((data.shape)) # first determine radii of all pixels
    r   = np.sqrt((x-center[0])**2+(y-center[1])**2)
    ind = np.argsort(r.flat) # get sorted indices
    sr  = r.flat[ind] # sorted radii
    sim = data.flat[ind] # image values sorted by radii
    ri  = sr.astype(np.int32) # integer part of radii (bin size = 1)
    csim = np.cumsum(sim, dtype=np.float64) # cumulative sum to figure out sums for each radii bin
    # determining distance between changes
    deltar = ri[1:] - ri[:-1] # assume all radii represented
    rind = np.where(deltar)[0] # location of changed radius

    if enclosed==1: #radially enclosed
        radialprofile = csim[rind[1:]]
    else: #radial averaging
        nr   = rind[1:] - rind[:-1] # number in radius bin
        tbin = csim[rind[1:]] - csim[rind[:-1]] # sum for image values in radius bins
        radialprofile = tbin/nr # the answer
    
    return radialprofile

def center(A, n):
    '''crop or embed A into nxn array from the center'''
    indx=(A.shape[0]-n+1)>>1
    indy=(A.shape[1]-n+1)>>1
    if indx >= 0 and indy >= 0:
        A2=A[indx:indx+n, indy:indy+n]
    elif indx <0 and indy <0 :
        A2=np.zeros((n,n))
        A2[-indx:-indx+A.shape[0], -indy:-indy+A.shape[1]]=A
    return A2

def photon_flux(magnitude, wvls):
    '''Claculate photon flux for magnitude at wavelength wvls'''
    Jy=1.51e7 #photons /sec / m^2 / dl_l
    name =  'UBVRIJHK';
    wvlc= np.array([0.35,0.44,0.55,0.64,0.79,1.26,1.60,2.22 ]) #center wavelength in micron
    dl_l= np.array([0.15,0.22,0.16,0.23,0.19,0.16,0.23,0.23 ])
    flux0=np.array([1810,4260,3640,3080,2550,1600,1080,670 ]) #zero magnitude flux
    flux= np.zeros((len(magnitude), len(wvls)))
    if np.max(wvls)<0.1:
        wvls=wvls*1e6
        
    for iwvl in range(len(wvls)):
        ind=np.argmin(abs(wvlc-wvls[iwvl]))
        zs=flux0[ind]*Jy*dl_l[ind]
        for imag in range(len(magnitude)):
            flux[imag, iwvl]=zs*10.**(-0.4*magnitude[imag]) #photon m^-2 s^-1
        #pde1=flux*dsa^2*dt*thru_tot*strehl;
    if len(flux)==1:
        flux=flux[0,0]
    return flux

def fftshift(a):
    a=a+0 #make a copy
    if a.dtype.name=='object':
        for it in np.ndindex(a.shape):
            a[it]=fftshift(a[it])
    else:
        na=a.ndim
        ia=na-2
        if ia<0:
            ia=0;
        for ic in range(ia, na):
            a=np.roll(a, int(a.shape[ic]/2), axis=ic)
    return a
def calc_psd(v, dt=1, fold=1, axis=-1):
    '''compute PSD from time histories'''
    nstep=v.shape[axis]
    nmid=(nstep>>1)
    df=1/(nstep*dt);
    af=np.abs(np.fft.fft(v,axis=axis))**2*(1/(df*nstep*nstep))
    
    if fold:
        f=np.arange(nmid+1)*df
        psd=af[:,0:1+nmid]
        psd[:,1:-1]+=af[:,-1:nmid:-1]
    else:
        f=np.arange(nstep)*df
        psd=af
    
    return(f,psd)
def plot_psd_cumu(f, psd, plot_options='-'):
    '''Plot the PSD cumulative integral from high to low'''
    fr=f[::-1]
    if psd.ndim==2:
        fr2=fr[:,None]
    else:
        fr2=fr;
    psdr=mysqrt(-cumtrapz(psd[::-1],fr2,axis=0,initial=0))
    loglog(fr,psdr, plot_options)

def cog(data):
    '''Center of gravity'''
    nx,ny=data.shape
    x=np.arange(nx)
    y=np.arange(ny)
    [Y,X]=np.meshgrid(y,x)
    ss=np.sum(data)
    cx=np.sum(X*data)/ss-(nx-1)/2
    cy=np.sum(Y*data)/ss-(ny-1)/2
    return cx,cy

def cog_shift(data):
    '''Shift cog to center'''
    cx,cy=cog(data)
    ccx=-int(round(cx))
    ccy=-int(round(cy))

    return np.roll(data, (ccx,ccy),axis=(0,1))
    
def plot_circle(radius, *args):
    '''plot_circle(radius, *args):'''
    theta=np.arange(0,np.pi*2,0.01)
    if type(radius) is list or type(radius) is np.ndarray:
        rx=radius[0]
        ry=radius[1]
    else:
        rx=radius
        ry=radius
    plot(rx*np.cos(theta), ry*np.sin(theta), *args)
    
def calc_width_gauss(dx,data):
    '''Compute the Gaussian width'''
    n=data.shape[0]
    n2=int(n/2)
    x=np.arange(-n2,n2)*dx
    XX,YY=np.meshgrid(x,x,indexing='xy')
    ds=np.sum(data)
    xb=np.sum(XX*data)/ds;
    yb=np.sum(YY*data)/ds;
    xsq=np.sum((XX-xb)**2*data)/ds
    ysq=np.sum((YY-yb)**2*data)/ds
    xysq=np.sum((XX-xb)*(YY-yb)*data)/ds
    #(ISO11146)
    mr=2*(xsq*ysq-xysq*xysq)**(1./4)
    return mr

def calc_fwhm(dx, intensity):
    '''A simple way to compute the FWHM'''
    return sqrt(np.sum(intensity>=0.5*np.max(intensity))*4/np.pi)*dx
