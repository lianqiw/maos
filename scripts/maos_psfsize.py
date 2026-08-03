import os
from pathlib import Path
import numpy as np
import glob


from readbin import readbin as read
#from libaos import read #this uses maos libaos.so
from astropy.modeling.models import Moffat2D, Const2D
from astropy.modeling.fitting import LevMarLSQFitter #handles bounded optim poorly
from astropy.modeling.fitting import TRFLSQFitter #handles bounded optim better
from astropy.modeling import Fittable2DModel, Parameter
from scipy.special import betaincinv
from scipy.integrate import dblquad
from scipy.optimize import brentq
from cycler import cycler

try:
    from natsort import natsorted
except:
    natsorted=sorted
import re

def rms(*args, **kargs):
    """compute rms of input"""
    if len(args) > 1:
        arr = np.array(args)
    else:
        arr = np.array(args[0])
    return np.sqrt(np.mean(np.real(arr * np.conj(arr)), **kargs))
def auto_crop_roi(img, threshold_rel=0.01, max_radius=None):
    ny, nx = img.shape
    if max_radius is None:
        max_radius=int(min(nx,ny)/2)
    # --- find peak
    iy, ix = np.unravel_index(np.argmax(img), img.shape)

    peak = img[iy, ix]
    thresh = peak * threshold_rel

    # --- radial growth (cheap approximation)
    r = 1
    while r < max_radius:
        y0, y1 = max(0, iy-r), min(ny, iy+r+1)
        x0, x1 = max(0, ix-r), min(nx, ix+r+1)

        sub = img[y0:y1, x0:x1]

        # stop when edges fall below threshold
        if np.all(sub[[0, -1], :] < thresh) and np.all(sub[:, [0, -1]] < thresh):
            break

        r += 1

    # final crop with margin
    margin = int(r * 0.3)
    y0 = max(0, iy - r - margin)
    y1 = min(ny, iy + r + margin + 1)
    x0 = max(0, ix - r - margin)
    x1 = min(nx, ix + r + margin + 1)

    return img[y0:y1, x0:x1], (x0, y0)

class EllipticalMoffat2D(Fittable2DModel):
    '''major-axis/axis-ratio parameterization of mofatt'''
    amplitude = Parameter(default=1.0)
    x_0 = Parameter(default=0.0)
    y_0 = Parameter(default=0.0)

    # Moffat scale along the major axis
    gamma = Parameter(default=1.0, bounds=(1e-6, None))

    # Minor / major axis ratio
    q = Parameter(default=1.0, bounds=(1e-3, 1.0))

    alpha = Parameter(default=2.5, bounds=(1e-6, None))

    # Rotation of major axis, radians
    theta = Parameter(default=0.0,
                      bounds=(-np.pi / 2, np.pi / 2))

    @staticmethod
    def evaluate(x, y, amplitude, x_0, y_0,
                 gamma, q, alpha, theta):

        dx = x - x_0
        dy = y - y_0

        cos_t = np.cos(theta)
        sin_t = np.sin(theta)

        # Rotate into the ellipse coordinate system
        xp = dx * cos_t + dy * sin_t
        yp = -dx * sin_t + dy * cos_t

        r2 = (xp / gamma)**2 + (yp / (q * gamma))**2

        return amplitude * (1.0 + r2)**(-alpha)
    
def fit_moffat_roi(img, alpha=3, name=None, circular=0):
    if True: #speed up
        sub, (xoff, yoff) = auto_crop_roi(img)
    else:
        sub=img
        xoff=0
        yoff=0
    ny, nx = sub.shape
    y, x = np.mgrid[:ny, :nx]

    # background
    #B0 = np.percentile(sub, 10)
    #B0=0
    # center (local)
    iy, ix = np.unravel_index(np.argmax(sub), sub.shape)
    if circular: #axial symmetric
        model = Moffat2D(
            amplitude=sub.max() ,#- B0,
            x_0=ix,
            y_0=iy,
            gamma=min(nx, ny) / 4,
            alpha=alpha
        ) #+ Const2D(amplitude=B0)
    else: #elliptical
        model = EllipticalMoffat2D(
        amplitude=np.max(sub),
        x_0=ix,
        y_0=iy,
        gamma=3.0,
        q=0.8,
        alpha=2.5,
        theta=0.0,
    )

    # bounded optimization
    model.alpha.bounds = (1+1e-6, 20) #below 1 not integratable, >>10 degenerate with gamma
    fitter = TRFLSQFitter()
        
    fit = fitter(model, x, y, sub)
    ierr= fitter.fit_info["status"]<1
    #Sanity check results
    model_img=fit(x,y)
    diff=np.sum((sub-model_img)**2)/np.sum(sub**2)
    if ierr or diff>0.05 or fit.alpha.value<=1+2e-7:
        print(f"Fitting is not good: {name}, diff={diff:.4g}, alpha={fit.alpha.value:.4g}, gamma={fit.gamma.value:.4g}, q={fit.q.value:.4g}")
        if ierr>4:
            print("nfev:", fitter.fit_info["nfev"])
            print("message:", fitter.fit_info["message"])
    
    # shift back to global coords
    fit.x_0.value += xoff
    fit.y_0.value += yoff
    
    return fit

        
def moffat_ensquare_width(gamma, alpha, f, q=1):
    def square_ee(W, gamma, alpha):
        half = W / 2
    
        def I(y, x):
            r2 = x*x + y*y
            return (1 + r2 / gamma**2) ** (-alpha)
    
        num, _ = dblquad(
            I,
            -half, half,
            lambda x: -half,
            lambda x: half
        )
    
        # total flux (analytic)
        total = np.pi * gamma**2 / (alpha - 1)
    
        return num / total

    def err(W):
        return square_ee(W, gamma, alpha) - f

    return brentq(err, 0.01*gamma, 50*gamma)*np.sqrt(q)
def moffat_slit_width(gamma, alpha, f, q=1):
    """
    Full slit width W enclosing fraction f for Astropy Moffat2D.

    Astropy convention:
        I(r) = A * [1 + (r/gamma)^2]^(-alpha)

    Parameters
    ----------
    gamma : float
        Core width parameter
    alpha : float
        Power-law index (must be > 1)
    f : float
        Desired enclosed energy fraction (0 < f < 1)

    Returns
    -------
    W : float
        Full slit width
    """
    if alpha <= 1:
        print(f"alpha={alpha}")
        raise ValueError("alpha must be > 1 for finite total energy")
    if not (0 < f < 1):
        raise ValueError("f must be between 0 and 1")

    # invert regularized incomplete beta
    u = betaincinv(0.5, alpha - 1.0, f)

    # convert to slit width
    a = gamma * np.sqrt(u / (1.0 - u))
    return 2.0 * a * np.sqrt(q)
def moffat_fwhm(gamma, alpha, q=1):
    return 2 * gamma * np.sqrt(2**(1/alpha)-1)*np.sqrt(q)
def moffat_encircle_width(gamma, alpha, frac, q=1):
    return 2 * gamma * np.sqrt((1 - frac)**(1/(1 - alpha)) - 1) * np.sqrt(q)
def parse_header_float(headers, key):
    res=[]
    for hh in headers:
        dp=''
        if isinstance(hh, str):
            s1=hh.find(f'{key} ')
            if s1!=-1:
                dp=hh[s1+10:]
            
        elif isinstance(hh, dict):
            dp=hh[key]
        if dp.find('/')==-1:
            print(dp)
            raise(Exception('Unable to parse header'))
        dp=float(dp[:dp.find('/')])
        res.append(dp)
    return np.stack(res)

def proc_psf(fn, fn_cache=None, **kargs):
    """
        Read PSF from a file and compute FWHM using Moffat fitting
    """
    if isinstance(fn, str):
        fn0=fn
        if fn_cache is None:
            fn_cache=fn0+'.npz' #cache results
    else:
        fn0=fn[0]
    #if fn_cache is None:
    #    print("Will not use cache for array if fn_cache is not set")
    skip=0
    if fn_cache and os.path.exists(fn_cache) and os.path.getmtime(fn_cache) > os.path.getmtime(fn0):
        try:
            data=np.load(fn_cache)
            dps=data["dps"]
            wvls=data["wvls"]
            alpha=data["alpha"]
            gamma=data["gamma"]
            if 'q' in data:
                q=data["q"]
            else:
                q=np.ones(gamma.size)
            #ress=data["ress"]
            skip=1
            #print(f'Using cached results from {fn_cache}')
        except:
            print(f'Loading cached results from {fn_cache} failed')
            pass

    if skip==0:
        datas=read(fn0)
        headers=[data.header for data in datas]
        if not isinstance(fn, str):
            count=1
            for fn2 in fn[1:]:
                datas+=read(fn2)
                count+=1
            if count>1:
                datas/=count
        dps=parse_header_float(headers, 'DP') #pixel width
        wvls=parse_header_float(headers, 'WVL')*1e6 #wavelength, convert to micron
        nwvl=datas.shape[0]
        gamma=np.zeros((nwvl))
        alpha=np.zeros((nwvl))
        q=np.zeros((nwvl))
        
        for iwvl in range(nwvl):
            #res=calc_fwhm_gaussian(datas[iwvl], dx=dps[iwvl])
            #res=maos_utils.print_psf_metrics(directory=f"{fd}/", x=-90, y=0, ee=50, seed=1)
            model=fit_moffat_roi(datas[iwvl], name=fn0, **kargs)
            alpha[iwvl]=model.alpha.value 
            gamma[iwvl]=model.gamma.value
            if "q" in model.param_names:
                q[iwvl]=model.q.value #elliptical
            else:
                q[iwvl]=1 #circular
        np.savez(fn_cache, dps=dps, wvls=wvls, alpha=alpha, gamma=gamma, q=q)
    #no longer caching ress
    nwvl=wvls.size
    ress=np.zeros((nwvl,3))
    for iwvl in range(nwvl):
        ress[iwvl, 0]=moffat_fwhm(gamma[iwvl], alpha[iwvl], q=q[iwvl])*dps[iwvl] #FWHM
        #ress[iwvl, 1]=moffat_encircle_width(gamma[iwvl], alpha[iwvl], 0.5)*dps[iwvl] #EE-50 diameter
        #ress[iwvl, 2]=moffat_encircle_width(gamma[iwvl], alpha[iwvl], 0.8)*dps[iwvl] #EE-80 en-circled diameter
        #ress[iwvl, 1]=moffat_ensquare_width(gamma[iwvl], alpha[iwvl], 0.5)*dps[iwvl] #EE-50 en-squared diameter
        #ress[iwvl, 2]=moffat_ensquare_width(gamma[iwvl], alpha[iwvl], 0.8)*dps[iwvl] #EE-80 en-squared diameter
        ress[iwvl, 1]=moffat_slit_width(gamma[iwvl], alpha[iwvl], 0.5, q=q[iwvl])*dps[iwvl] #EE-50 slit width
        ress[iwvl, 2]=moffat_slit_width(gamma[iwvl], alpha[iwvl], 0.8, q=q[iwvl])*dps[iwvl] #EE-80 slit width
        #alpha is manually set to 3, so the ratio between fwhm and EE-p% is constant
        #not a good idea to fit over alpha
    
    #print(f'Save results to {fn_cache}')
    return ress,dps,wvls

def proc_psfs(fd, seeds=[1], **kargs):
    """
    Read PSFs from a folder and compute FWHM.
    Multiple seeds are averaged before fitting
    """
    ressc=[] 

    ress=[] #result array
    fr=[] #field point radius
    fx=[] #field point x
    fy=[] #field point y
    dps=None #PSF pixel sampling
    fns=natsorted(glob.glob(f"{fd}/evlpsfcl_{seeds[0]}_x*_y*.fits"))
    for fn in fns:
        m = re.search(r"_x([+-]?\d+(?:\.\d+)?)_y([+-]?\d+(?:\.\d+)?)\.fits$", fn)
        x, y=m.groups()
        fns2=[fn]
        for seed in seeds[1:]:
            fns2.append(f"{fd}/evlpsfcl_{seed}_x{x}_y{y}.fits")
        x, y = map(float, m.groups())
        fr.append(np.sqrt(x*x+y*y))
        fx.append(x)
        fy.append(y)
        fn_cache=f"{fd}/evlpsfcl_{'_'.join(map(str, seeds))}_x{x}_y{y}.fits.npz"
        res,dps,wvls=proc_psf(fns2, fn_cache, **kargs)
        ress.append(res)
    fr=np.array(fr)
    fx=np.array(fx)
    fy=np.array(fy)

    ress=np.stack(ress)
    ind=np.lexsort((fy, fx))
    fr=fr[ind]
    fx=fx[ind]
    fy=fy[ind]
    ress=ress[ind]
    avg=np.mean(ress,axis=0) #field averaged result
    err=np.std(ress,axis=0) #variation over the field
    
    #open loop
    fn=f"{fd}/evlpsfol_{seeds[0]}.fits"
    fns2=[fn]
    for seed in seeds[1:]:
        fns2.append(f"{fd}/evlpsfol_{seed}.fits")
    fn_cache=f"{fd}/evlpsfol_{'_'.join(map(str, seeds))}.fits.npz"
    ressol, dps, wvls=proc_psf(fns2, fn_cache, **kargs)

    ressc.append({'avg':avg, 'std':err, 'cl':ress, 'ol':ressol,'fr':fr, 'fx':fx, 'fy':fy, 'dps':dps, 'wvls':wvls})
    if len(ressc)>1: #collect seeds
        ressc = {k: rms(np.array([d[k] for d in ressc if d is not None]), axis=0) for k in ressc[0]}
    elif len(ressc)==1:
        ressc = ressc[0]
    else:
        ressc = None
    return ressc
    
