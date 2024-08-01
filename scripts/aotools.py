#!/usr/bin/env python3

#Helper routines to parse MAOS results and data processing

import glob
import os
import numpy as np
import scipy
import warnings
np.set_printoptions(threshold=100,suppress=False,precision=4,floatmode='maxprec',linewidth=120)
from scipy.special import erf
from scipy.integrate import cumtrapz
from numpy import sqrt, exp, log, floor, ceil, nan, pi, sin, cos, tan, arcsin, arccos
from numpy.random import rand, randn
from numpy.fft import fft,ifft,fft2, ifft2, fftshift
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, plot, semilogx, semilogy, loglog, xlabel, ylabel, legend, grid, clf, subplot, xlabel, ylabel, title, xlim, ylim,clim, close
from cycler import cycler
try:
    from natsort import natsorted
except:
    natsorted=sorted
try:
    import libaos as aos
except:
    import lib2py #it generates libaos
    import libaos as aos
from libaos import read,write
from readbin import readbin, get_header, header
from draw import draw, locembed, coord2grid
from maos_result import maos_res, maos_res_tot, maos_res_hi, maos_res_each, maos_cumu
#import maos_client
def error(msg):
    raise(Exception(msg))
def sind(x):
    return np.sin(x*np.pi/180)

def cosd(x):
    return np.cos(x*np.pi/180)

def tand(x):
    return np.tan(x*np.pi/180)

def iscell(arr):
    return type(arr)==np.ndarray and arr.dtype.name=='object'

'''use of .flat with an FORTRAN ordered array will lead to non-optimal memory
access as adjacent elements in the flattened array (iterator, actually) are not
contiguous in memory'''

def isequal(a, b):
    if type(a)==np.ndarray and type(b)==np.ndarray:
        if iscell(a):
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

def mysqrt(x):
    if type(x) is np.ndarray:
        return np.sign(x)*np.sqrt(np.abs(x))
    elif x<0:
        return -np.sqrt(-x)
    else:
        return np.sqrt(x)

def rss(*args, **kargs):
    '''compute rss of input'''
    if len(args)>1:
        arr=np.array(args)
    else:
        arr=np.array(args[0])
    return mysqrt(np.sum(np.real(np.sign(arr)*arr*np.conj(arr)), **kargs))
def rms(*args, **kargs):
    '''compute rms of input'''
    if len(args)>1:
        arr=np.array(args)
    else:
        arr=np.array(args[0])
    return np.sqrt(np.mean(np.real(arr*np.conj(arr)), **kargs))

def styles(ii):
    reset_color()
    lines=['-','--','-.',':'];
    return lines[np.mod(ii,len(lines))]

def reset_color():
    plt.gca().set_prop_cycle(None)

def cycler_marker(line='-'):
    markercycle = cycler(marker='+x*.Xodx*.Xod+') #['+', 'x', '*', '.', 'X','o','d'])
    colorcycle = cycler(color='bgrcmykbgrcmyk')
    linecycle= cycler(linestyle=[line])
    plt.gca().set_prop_cycle((colorcycle + markercycle)*linecycle) # gca()=current axis

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
#def figure(*args, **kargs):
#    plt.figure(*args, **kargs)
#    plt.clf()


def width_at(x, y, height):
    '''compute width of 1d function defined on x above a threshold height'''
    return len(np.argwhere(y>=max(y)*height))*(x[1]-x[0])

def cellsum(x):
    xsum=np.zeros(x.shape)
    for ii in range(x.size):
        jj=np.unravel_index(ii, x.shape)
        xsum[jj]=x[jj].sum()

    return xsum

def cellmm(A,B):
    '''
        Multiply matrices A and B using MAOS convention
        A or B may be cell or dense matrices.
        Python stores matrix in row major order, which is like everything must be transposed.
    '''
    while A.ndim>2 and A.shape[0]==1:
        A=A[0]
    while B.ndim>2 and B.shape[0]==1:
        B=B[0]
    
    if A.ndim>2 or B.ndim>2:
        error("Invalid input for this operation")
    if A.ndim<2:
        A=A.reshape((1, A.size))
    if B.ndim<2:
        B=B.reshape((1, B.size))
    if iscell(A) and iscell(B):
        C=np.empty((B.shape[0], A.shape[1]), dtype=object)
        for iz in range(C.shape[0]):
            for ix in range(C.shape[1]):
                for iy in range(A.shape[0]):
                    if A[iy,ix] is None or A[iy,ix].size==0 or B[iz,iy] is None or B[iz,iy].size==0:
                        continue
                    print(ix,iy,iz,A[iy,ix].shape, B[iz,iy].shape)
                    if C[iz,ix] is None:
                        C[iz,ix]=B[iz,iy]@A[iy,ix]
                    else:
                        C[iz,ix]+=B[iz,iy]@A[iy,ix]
                    #print(ix,iy,iz,A[iy,ix].shape, B[iz,iy].shape, C[iz,ix].shape)
    else:
        C=B@A
    return C
def celltrans(A):
    '''Transpose cell arrays'''
    if A.ndim==1:
        A=np.reshape(A, (1, A.size))
    AT=np.copy(A.T)
    if iscell(A):
        for i in range(AT.size):
            AT.flat[i]=np.transpose(AT.flat[i])
    return AT

#remove tip/tilt/focus from gradients
def grad_ttfr(grad, saloc):
    '''remove tip/tilt/focus move from gradients defined on subaperture location saloc'''
    if iscell(grad):
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
def loc_zernike_proj(loc, rmin, rmax, radonly=0):
    '''Project onto zernike mode between rmin and rmax from opd which is defined on loc'''
    #D=np.max(np.max(loc,axis=1)-np.min(loc, axis=1))
    mod=aos.zernike(loc, 0, rmin, rmax, radonly).T
    rmod=np.linalg.pinv(mod)
    return mod, rmod
#remove zernike modes from rmin to rmax from 1-D OPD and loc
def opd_loc_project_zernike(opd, loc, rmin, rmax, radonly=0):
    '''Project onto zernike mode between rmin and rmax from opd which is defined on loc'''
    mod, rmod=loc_zernike_proj(loc, rmin, rmax, radonly)
    return mod@(rmod@opd)

def opd_loc_remove_zernike(opd, loc, rmin, rmax, radonly=0):
    '''Remove zernike mode between rmin and rmax from opd which is defined on loc'''
    mod, rmod=loc_zernike_proj(loc, rmin, rmax, radonly)
    opd2=opd-mod@(rmod@opd)
    return opd2
def mk2dloc(shape):
    nx=shape[0]
    if len(shape)>1:
        ny=shape[1]
    else:
        ny=nx
    return aos.mksqloc(ny, nx, 1, 1, -ny/2, -nx/2).reshape(2,nx*ny)
def mktt(amp):
    '''create piston/tip/tilt mode from amplitude map'''
    (nx,ny)=amp.shape
    loc=aos.mksqloc(ny, nx, 2./ny, 2./nx, -1, -1).reshape(2, amp.size)
    loc[:,amp.flat<=0]=0
    return loc
def mkptt(amp):
    '''create piston/tip/tilt mode from amplitude map'''
    (nx,ny)=amp.shape
    loc=aos.mksqloc(ny, nx, 1, 1, -ny/2, -nx/2).reshape(2, amp.size)
    pis=(amp>0).reshape(1,amp.size)
    loc[:,pis[0]==0]=0
    return np.r_[pis, loc]

#remove zernike modes from rmin to rmax from 2-D OPD
def opd_remove_zernike(opd, mask, rmin, rmax, radonly=0):
    '''OPD is 2d array, mask indicates valid points, rmin and rmax are minimum and maximum zernike order (inclusive)'''
    opdshape=opd.shape
    #oloc=mksqloc(opd.shape[1], opd.shape[0], 1, 1, -opd.shape[1]/2, -opd.shape[0]/2).reshape(2,opd.size)
    oloc=mk2dloc(opd.shape)
    opd=opd.reshape((opd.size,))
    if mask is None:
        mask=opd!=0
    else:
        mask=mask.reshape((mask.size,))
    oloc2=oloc[:,mask].copy() #copy is necessary for C to access the data
    opd2=np.zeros(opd.shape)
    opd2[mask]=opd_loc_remove_zernike(opd[mask], oloc2, rmin, rmax, radonly);
    #D=np.max(np.max(oloc2,axis=1)-np.min(oloc2, axis=1))
    #mod=zernike(oloc2, D, rmin, rmax, radonly).T
    #rmod=np.linalg.pinv(mod)
    #opd2[mask]=opd[mask]-mod@(rmod@opd[mask])
    return opd2.reshape(opdshape)
def opd_loc_remove_focus(opd, loc):
    '''Remove focus mode from opd which is defined on loc'''
    return opd_loc_remove_zernike(2,2,1)
def read_many(fdin):
    '''read many files together'''
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
    try:
        res=simplify(np.array(res, dtype=object))
    except:
        pass
    return res,fds
def read_many_dict(fdin):
    '''read many files together into a dictionary'''
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
    """create an identity matrix"""
    dd=np.zeros((nx,nx))
    np.fill_diagonal(dd, val)
    return dd
def svd_inv(A, thres=0, tikcr=0):
    '''SVD inversion of matrix: svd_inv(A, thres=0, tikcr=0)'''
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

def plot_cdf(y, *args, **kargs):
    '''plot cumulative distribution function. Take absolute value times 2 if abs is set'''
    n=y.size
    x=np.arange(n)/(n-1)
    if 'abs' in kargs and kargs['abs']:
        y=np.abs(y)*2
    plot(np.sort(y.flat), x, *args, **kargs)
    ylim(0,1)

def plot_cross(A, n=0, dx=1):
    '''Plot a cross section'''
    if n==0 or n is None:
        n=A.shape[0]
        
    n1=max((A.shape[0]-n)>>1,0)
    n2=min((A.shape[0]+n)>>1,A.shape[0])
    print(n1,n2,A.shape[0]>>1)
    plot((np.arange(n1,n2)-(A.shape[0]>>1))*dx, A[n1:n2,A.shape[1]>>1])

def plot_smooth(x,y,*args,**kargs):
    from scipy.interpolate import make_interp_spline, BSpline
    ind=np.argsort(x)
    x=x[ind]
    y=y[ind]
    #define x as 200 equally spaced values between the min and max of original x
    xnew = np.linspace(x.min(), x.max(), 200)
    #define spline
    spl = make_interp_spline(x, y, k=3, bc_type="natural")
    y_smooth = spl(xnew)

    #create smooth line chart
    plt.plot(xnew, y_smooth, *args, **kargs)
def radial_profile(data, nx=None):
    if nx is None:
        nx=data.shape[0]>>1
    xx=np.arange(nx).astype(np.float64).copy() # np.linspace(0,nx-1,nx) #np.arange(nx).asfarray()
    return aos.calcenc(data, xx, -1, 10) #azimuthal averaging
def radial_profile_obsolete(data, center=None, enclosed=0):
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
def check2d(A, nx=None): 
    '''check shape of A and return a 2-d square array'''
    if A is None:
        return A
    while iscell(A) or A.ndim>2:
        if A.shape[0]==1:
            A=A[0]
        else:
            error('Unable to process A with shape ', A.shape)
    if A.ndim==2 and A.shape[0]==A.shape[1]:
        return A
    if nx is None:
        nx=round(np.sqrt(A.size))
    ny=round(A.size/nx)
    if A.size!=nx*ny:
        error('Please provide correct nx')
    return A.reshape(nx,ny)

def center(A, nx=None, ny=None, ratio=None):
    if type(A) is list:
        B=[]
        for Ai in A:
            B.append(center(Ai, nx, ny, ratio))
        return B
    elif A is None:
        return None
    '''crop or embed A into nxn array from the center'''
    A=check2d(A)
    if ratio is not None:
        nx=int(A.shape[0]*ratio)
        ny=int(A.shape[1]*ratio)
    else:
        nx=int(nx)
        if ny is None:
            ny=int(np.round(A.shape[1]*nx/A.shape[0]))
        else:
            ny=int(ny)
    if A.shape[0]>=nx and A.shape[1]>=ny:#extract
        indx=(A.shape[0]-nx+1)>>1
        indy=(A.shape[1]-ny+1)>>1
        A2=A[indx:indx+nx, indy:indy+ny]
    elif A.shape[0]<nx and A.shape[1]<ny:#embed
        indx=(nx-A.shape[0]+1)>>1
        indy=(nx-A.shape[1]+1)>>1
        A2=np.zeros((nx,ny),dtype=A.dtype)
        A2[indx:indx+A.shape[0], indy:indy+A.shape[1]]=A
    else:
        raise(Exception('Does not support this operation'))
    return A2

def embed(a,ratio=2):
    '''Embed a 2-d image into a new array bigger by ratio '''
    return center(a, ratio=ratio)

def unembed(a,ratio=2):
    '''undo embed(). unembed center piece'''
    return center(a, ratio=1/ratio)

def photon_flux(magnitude, wvls):
    '''Claculate photon flux for magnitude at wavelength wvls'''
    Jy=1.51e7 #photons /sec / m^2 / dl_l
    name =  'UBVRIJHK'
    wvlc= np.array([0.36,0.44,0.55,0.64,0.79,1.26,1.60,2.22 ]) #center wavelength in micron
    dl_l= np.array([0.15,0.22,0.16,0.23,0.19,0.16,0.23,0.23 ])
    flux0=np.array([1810,4260,3640,3080,2550,1600,1080,670 ]) #zero magnitude flux
    magnitude=np.asarray(magnitude)
    magnitude.shape=(magnitude.size,)
    wvls=np.asarray(wvls)
    wvls.shape=(wvls.size,)
    nmag=magnitude.size
    nwvl=wvls.size
    flux= np.zeros((nmag, nwvl))
    if np.max(wvls)<0.1:
        wvls=wvls*1e6

    for iwvl in range(nwvl):
        ind=np.argmin(abs(wvlc-wvls[iwvl]))
        zs=flux0[ind]*Jy*dl_l[ind]
        for imag in range(nmag):
            flux[imag, iwvl]=zs*10.**(-0.4*magnitude[imag]) #photon m^-2 s^-1
        #pde1=flux*dsa^2*dt*thru_tot*strehl;
    if flux.size==1:
        flux=flux[0,0]
    return flux

def calc_psd(v, dt=1, fold=1, axis=-1):
    '''compute PSD from time histories'''
    nstep=v.shape[axis]
    nmid=(nstep>>1)
    df=1/(nstep*dt);
    af=np.abs(fft(v,axis=axis))**2*(1/(df*nstep*nstep))

    if fold:
        f=np.arange(nmid+1)*df
        if af.ndim==2:
            af[:,1:nmid]+=af[:,-1:nmid:-1]
            psd=af[:,0:1+nmid] #psd[:,1:-1]+=af[:,-1:nmid:-1]
        elif af.ndim==1:
            af[1:nmid]+=af[-1:nmid:-1]
            psd=af[0:1+nmid] #psd[:,1:-1]+=af[:,-1:nmid:-1]
        else:
            raise(Exception('Please fill in'))
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

def cog(data, center=None):
    '''Center of gravity'''
    nx,ny=data.shape
    x=np.arange(nx)
    y=np.arange(ny)
    if center is None:
        center=((nx-1)/2,(ny-1)/2)
    [Y,X]=np.meshgrid(y,x)
    ss=np.sum(data)
    cx=np.sum(X*data)/ss-center[0]
    cy=np.sum(Y*data)/ss-center[1]
    return cx,cy

def cog_shift(data, center=None):
    '''Shift cog to center'''
    cx,cy=cog(data, center)
    ccx=-int(round(cx))
    ccy=-int(round(cy))

    return np.roll(data, (ccx,ccy),axis=(0,1))

def plot_circle(radius, center=[0, 0], *args):
    '''plot_circle(radius, *args):'''
    theta=np.arange(0,np.pi*2,0.01)
    if isinstance(radius, (list, tuple, np.ndarray)):
        rx=radius[0]
        ry=radius[1]
    else:
        rx=radius
        ry=radius
    if isinstance(center, (list, tuple, np.ndarray))==False:
        center=[center,center]
    plot(rx*np.cos(theta)+center[0], ry*np.sin(theta)+center[1], *args)

def calc_width_gauss(dx,data):
    '''Compute the Gaussian width of 2-d array data with sampling of dx'''
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


def stfun(opd, amp):
    '''Compute structure function.'''
    nx=opd.shape[0]
    if amp is None:
        amp=np.ones(opd.shape)
    opd2=center(opd, nx*2);
    amp2=center(amp, nx*2);
    hat0=fft2(amp2);
    hat1=fft2(opd2*amp2);
    hat2=fft2(opd2*opd2*amp2);
    hat3=2*(np.real(hat2*np.conj(hat0))-hat1*np.conj(hat1))
    hat3i=np.real(ifft2(hat3))
    hat0i=np.real(ifft2(hat0*np.conj(hat0)))
    hat3i[hat0i>0.4]/=hat0i[hat0i>0.4]
    hat3i[hat0i<=0.5]=0;
    st=fftshift(hat3i);
    return st
def opdfillzero(im0, nrep):
    '''Fill zero value of 2d image im0 nrep times by shifting and adding'''
    im=im0+0 #avoid changing input data
    for irep in range(nrep):
        mask=im==0;#select zero points to fill
        print(np.sum(mask))
        if np.sum(mask)==0:
            break;

        mask1=im!=0;
        im2=0;
        mask2=0;
        for ix in [-1,1]:
            for axis in 0,1:
                im2+=np.roll(im, ix, axis=axis)
                mask2+=np.roll(mask1, ix, axis=axis)
        mask[mask2==0]=0
        im[mask]=im2[mask]/mask2[mask]
    return im

def str_remove_prefix(ini_strlist, offset=0):
    '''find and remove the common prefix of string list'''
    if len(ini_strlist)>1:
        '''find the common prefix of string list'''
        prefix = ini_strlist[0]
        for string in ini_strlist[1:]:
            while string[:len(prefix)] != prefix and prefix:
                prefix = prefix[:len(prefix)-1]
            if not prefix:
                break
        nprefix=len(prefix)+offset
        return [ tmp[nprefix:] for tmp in ini_strlist]
    else:
        return ini_strlist

def meshgrid(nx,dx):
    xv=(np.arange(nx)-nx/2)*dx
    [X,Y]=np.meshgrid(xv, xv)
    return X,Y

def reflection(n1, n2, degree):
    '''incident from n1 to n2 (refractive index), incident angle in degree
        compute s and p wave reflectivity
    '''
    ti=degree*np.pi/180 #incident angle in radian
    tt=np.arcsin(np.sin(ti)*n1/n2) #refractive angle in radian
    ci=np.cos(ti)
    ct=np.cos(tt)
    rs=((n1*ci-n2*ct)/(n1*ci+n2*ct))**2
    rp=((n1*ct-n2*ci)/(n1*ct+n2*ci))**2
    return rs,rp

def make_R2(n, dx=1):
    '''Create radius squared of a meshgrid'''
    X=np.arange(-n/2,n/2)*dx
    X.shape=(n,1)
    Y=X+0
    Y.shape=(1,n)
    R2=X*X+Y*Y
    return R2
def gen_atm(nx, dx=1./64., r0=0.186, L0=30., powerlaw=-11./3., seed=0):
    '''Generate atmosphere'''
    R2=make_R2(nx) #radius squared (with dx=1)
    if L0>0: #outer scale
        R2+=(nx*dx/L0)**2
    if seed:
        rng=np.random.default_rng(seed)
    else:
        rng=np.random.default_rng()
    wn=rng.standard_normal(nx*nx)+1j*rng.standard_normal(nx*nx)
    wn.shape=(nx,nx)
    alpha=0.1517*(nx*dx/r0)**(5./6.)*0.5e-6/(2*np.pi)
    spect=alpha*wn*R2**(-abs(powerlaw/4)) #spectrum
    spect[nx>>1, nx>>1]=0
    spect[:,0]=0
    spect[0,:]=0
    screen=np.real(fft2(fftshift(spect)))
    return screen
def abs2(A):
    return np.real(A*np.conj(A))

def myfft2(A):
    return fftshift(fft2(fftshift(A)))
def myifft2(A):
    return fftshift(ifft2(fftshift(A)))
    
def calc_psf(amp, opd, wvl, fembed=2, wts=None):
    '''
        Compute PSF normalized by Strehl ratio.
        For multi-wavelength PSF, the sampling is matched to the first wavelength
    '''
    
    amp=check2d(amp)
    opd=check2d(opd)
    if amp.ndim!=2 or amp.shape[0]!=amp.shape[1]:
        raise(Exception('Only supports square array'))
    if type(wvl) is list:
        wvl=np.asarray(wvl)
    if type(wvl) is np.ndarray:
        psfs=0
        wt=1/wvl.size

        nx=amp.shape[0]
        wvl0=np.mean(wvl) #.flat[0]
        if type(wts) is list:
            wts=np.asarray(wts)
        for iwvl in range(wvl.size):
            wvli=wvl.flat[iwvl]
            #nx2=int(round(wvli/wvl0*nx/4))*4 #scale array size so that the sampling match the first wvl.
            if type(wts) is np.ndarray:
                wt=wts.flat[iwvl]
            #psfs+=center(calc_psf(center(amp,nx2),center(opd,nx2),wvli),nx)*wt
            psfi=calc_psf(amp,opd,wvli,fembed=fembed*(wvli/wvl0))*wt
            #if psfs is None:
            #    psfs=psfi
            #else:
            psfs+=center(psfi,nx*fembed)
        return psfs
    else:
        if opd is None or opd.size==0:
            wvf=amp
        else:
            wvf=amp*np.exp(opd*(2j*pi/wvl))

        nx2=round(fembed*wvf.shape[0]/2)*2
        psf=abs2(myfft2(center(wvf,nx2)))*(1/np.sum(amp)**2)
        return psf

def gerchberg_saxton(A1, A2, P1B=None, modes=None, rmodes=None, fembed=2, nrep=100, doprint=0):
    '''
    The Gerchberg Saxton algorithm (1972)
    Parameters
    ------------
    A1: pupil plane amplitude, embeded into a square array to match PSF sampling
    A2: image plane amplitude (sqrt of PSF), embeded into a square array with the same size as AMP
    modes: if set, the estimated OPD will be projected onto these modes.
    P1B: the initial phase (in radian). If not set, will use random values. When projection onto the modes, the content is subtracted first.
    fembed: embedding factor ratio
    nrep: number of iterations
    doprint: plot the convergence and result
    '''
    err=np.zeros((nrep,))
    #rng = np.random.default_rng(1)
    if fembed!=1:
        A1=embed(A1, ratio=fembed)
    if A1.shape[0]!=A2.shape[0]: 
        A2=center(A2, A1.shape[0])
    if P1B is not None:
        P1B=center(P1B, A1.shape[0])

    if P1B is not None:
        P1E=P1B.copy()
    else:
        P1E=np.zeros(A1.shape)
        #P1E=(rng.random(A1.shape)-0.5)*np.pi/2 #initialize guess
        
    C1=A1*exp(1j*P1E)
    if modes is not None and rmodes is None: #projection solution to the modes
        rmodes=np.linalg.pinv(modes)
    alpha=1./np.sum(A2**2)
    A2=fftshift(A2)
    for i in range(nrep):
        C2=fft2(C1) #no need to apply normalization 
        if doprint:
            err[i]=alpha*np.sum((np.abs(C2)-A2)**2)
        C2*=A2/np.abs(C2) #fix amplitude
        C1=ifft2(C2) #no need to apply normalization 
        if modes is not None or i+1==nrep:
            P1E=np.angle(C1)
            P1E[A1>0]-=np.mean(P1E[A1>0]) #remove piston
            P1E[A1<=0]=0 #remove points out side of aperture
            
        if modes is None:
            C1*=A1/np.abs(C1) #fix amplitude
        else:
            if P1B is not None:
                P1E-=P1B
            P1E=unembed(P1E, ratio=fembed)
            mv=(P1E.flat@rmodes) #project onto modes
            P1E=np.reshape(mv@modes, P1E.shape) 
            P1E=embed(P1E, ratio=fembed)
            if P1B is not None:
                P1E+=P1B
            C1=A1*exp(1j*P1E)

    if doprint:
        print(err)
    if fembed!=1:
        P1E=unembed(P1E, ratio=fembed)
    if modes is None:
        return P1E
    else:
        return P1E,mv

def split_even_odd(A, doshift=0):
    '''
        split A into even and odd components. A has dc term at lower left corner A[0,0].
    '''
    if doshift:
        Ae=fftshift(A)
    else:
        Ae=A+0
    Ae[1:,1:]+=Ae[-1:0:-1,-1:0:-1]
    #Ae[:,0]=0; Ae[0,:]=0
    Ae*=0.5;
    if doshift:
        Ae=fftshift(Ae)
    Ao=A-Ae
    return Ae, Ao

def fast_furious(AMP, PSF, PSFD, OPDD, epsilon=0.1, doplot=0):
    '''
    The fast and furious algorithm for weak phase
    Parameters
    ------------
    AMP: pupil amplitude, embeded into a square array to match PSF sampling. max is 1.
    PSF: measured PSF with unknown OPD, embeded into a square array with the same size as AMP
    PSFD: measured PSF with OPD+OPDD, embeded into a square array with the same size as AMP
    OPDD: opd of phase diversity, embeded into a square array with the same size as AMP
    epsilon=athres: regularization     
    doplot: do plotting
    '''
    a=np.real(myfft2(AMP)) #do not use sqrt(PSF_DL) as a might be negative
    PSF_DL=a**2 #Diffraction limited PSF
    pp=PSF+(1-np.max(PSF)/np.max(PSF_DL))*PSF_DL #modified measurements with peak that matches PSF_DL 
    ppe,ppo=split_even_odd(pp) #split to even/odd components
    pd=PSFD+(1-np.max(PSFD)/np.max(PSF_DL))*PSF_DL #modified measurements with diversity
    pde,pdo=split_even_odd(pd)
    yp=-a*ppo/(2*PSF_DL+epsilon) #paper was wrong. missing - sign
    AMP_OPDDe, AMP_OPDDo=split_even_odd(OPDD*AMP)
    vd=np.real(myfft2(AMP_OPDDe))
    yd=np.imag(myfft2(AMP_OPDDo))
    #figure(); draw([vd, yd])
    vs=(pde-ppe-vd**2-yd**2-2*yp*yd)/(2*vd+np.max(vd)*epsilon) #paper was wrong in the ordering of the first two tems
    vp=np.sign(vs)*sqrt(np.abs(ppe-PSF_DL-yp**2))
    #figure();draw([vs,vp])
    
    #estimate
    athres=0.5
    P1o=np.real(myifft2(1j*yp))
    P1o[AMP>athres]/=AMP[AMP>athres]; P1o[AMP<=athres]=0; P1o[AMP>athres]-=np.mean(P1o[AMP>athres])
    P1e=np.real(myifft2(vp))
    P1e[AMP>athres]/=AMP[AMP>athres]; P1e[AMP<=athres]=0; P1e[AMP>athres]-=np.mean(P1e[AMP>athres])
    if doplot:
        figure(); draw([P1e, P1o, AMP_OPDDe, AMP_OPDDo])

    return P1o+P1e

def nanmean(A, *args, **kargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(A, *args, **kargs)
    
def nanmedian(A, *args, **kargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmedian(A, *args, **kargs)
    
def nanmin(A, *args, **kargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmin(A, *args, **kargs)
    
def dtf_blur(psf, notf, blur):
    '''Add (blur>0) or remove (blur<0) detector blury function'''
    if psf.ndim==3:
        psfout=np.zeros(psf.shape)
        for i in range(psf.shape[0]):
            psfout[i]=dtf_blur(psf[i], notf, blur)
        return psfout
    wvf=fft2(center(check2d(psf), notf))
    otf=aos.dtf_otf(notf, notf, 1, 1, 0, abs(blur), 2)
    #aos.deblur(psf2, otf)
    if blur>=0:
        wvf*=otf
    else:
        wvf[otf!=0]/=otf[otf!=0]
    psfout=np.real(ifft2(wvf))
    psfout[psfout<0]=0
    return center(psfout, psf.shape[0])


def interp_cn2(ht2, ht, wt):
    '''Represent the cn2dh in new layers by interpolating the forward and backward cumulative sum'''
    wt2=np.interp(ht2, ht, np.cumsum(wt))
    wt2rev=np.interp(ht2, ht, np.cumsum(wt[::-1])[::-1])

    wt2[1:]-=wt2[0:-1] #undo cumulative
    wt2rev[0:-1]-=wt2rev[1:] #undo reverse cumulative
    wt2=0.5*(wt2+wt2rev) #average
    return wt2/np.sum(wt2)

def strehl(rms, wvl):
    return exp(-(2*np.pi*rms/wvl)**2)

def savefig(fn, *args, **kargs):
    plt.gcf().tight_layout() #fix layout
    plt.savefig(fn, *args, **kargs, bbox_inches='tight')
