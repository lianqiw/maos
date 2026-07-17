import os
from pathlib import Path
import numpy as np
import glob


from readbin import readbin as read
#from libaos import read #this uses maos libaos.so
from astropy.modeling.models import Moffat2D, Const2D
from astropy.modeling.fitting import LevMarLSQFitter
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
def auto_crop_roi(img, threshold_rel=0.2, max_radius=50):
    ny, nx = img.shape

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
    
def fit_moffat_roi(img, alpha=3):
    sub, (xoff, yoff) = auto_crop_roi(img)

    ny, nx = sub.shape
    y, x = np.mgrid[:ny, :nx]

    # background
    #B0 = np.percentile(sub, 10)
    #B0=0
    # center (local)
    iy, ix = np.unravel_index(np.argmax(sub), sub.shape)

    model = Moffat2D(
        amplitude=sub.max() ,#- B0,
        x_0=ix,
        y_0=iy,
        gamma=min(nx, ny) / 4,
        alpha=alpha
    ) #+ Const2D(amplitude=B0)

    # stabilize
    model.alpha.fixed = True #for stability
    model.gamma.bounds = (1e-2, None)

    fitter = LevMarLSQFitter()
    fit = fitter(model, x, y, sub)

    # shift back to global coords
    fit.x_0.value += xoff
    fit.y_0.value += yoff

    return fit

        
def moffat_ensquare_width(gamma, alpha, f):
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

    return brentq(err, 0.01*gamma, 50*gamma)
def moffat_slit_width(gamma, alpha, f):
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
        raise ValueError("alpha must be > 1 for finite total energy")
    if not (0 < f < 1):
        raise ValueError("f must be between 0 and 1")

    # invert regularized incomplete beta
    u = betaincinv(0.5, alpha - 1.0, f)

    # convert to slit width
    a = gamma * np.sqrt(u / (1.0 - u))
    return 2.0 * a
def moffat_fwhm(gamma, alpha):
    return 2 * gamma * np.sqrt(2**(1/alpha)-1)
def moffat_encircle_width(gamma, alpha, frac):
    return 2 * gamma * np.sqrt((1 - frac)**(1/(1 - alpha)) - 1)
def parse_header_float(psfs, key):
    res=[]
    for psfi in psfs:
        hh=psfi.header
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

def proc_psf(fn, **kargs):
    """
        Read PSF from a file and compute FWHM using Moffat fitting
    """
    fn_cache=fn+'.npz' #cache results
    skip=0
    if os.path.exists(fn_cache) and os.path.getmtime(fn_cache) > os.path.getmtime(fn):
        try:
            data=np.load(fn_cache)
            dps=data["dps"]
            wvls=data["wvls"]
            alpha=data["alpha"]
            gamma=data["gamma"]
            ress=data["ress"]
            skip=1
        except:
            pass
    if skip==0:
        datas=read(fn)
        dps=parse_header_float(datas, 'DP') #pixel width
        wvls=parse_header_float(datas, 'WVL')*1e6 #wavelength, convert to micron
        nwvl=datas.shape[0]
        gamma=np.zeros((nwvl))
        alpha=np.zeros((nwvl))
        ress=np.zeros((nwvl,3))
        for iwvl in range(nwvl):
            #res=calc_fwhm_gaussian(datas[iwvl], dx=dps[iwvl])
            #res=maos_utils.print_psf_metrics(directory=f"{fd}/", x=-90, y=0, ee=50, seed=1)
            model=fit_moffat_roi(datas[iwvl], **kargs)
            alpha[iwvl]=model.alpha.value 
            gamma[iwvl]=model.gamma.value
            ress[iwvl, 0]=moffat_fwhm(gamma[iwvl], alpha[iwvl])*dps[iwvl] #FWHM
            #ress[iwvl, 1]=moffat_encircle_width(gamma[iwvl], alpha[iwvl], 0.5)*dps[iwvl] #EE-50 diameter
            #ress[iwvl, 2]=moffat_encircle_width(gamma[iwvl], alpha[iwvl], 0.8)*dps[iwvl] #EE-80 en-circled diameter
            #ress[iwvl, 1]=moffat_ensquare_width(gamma[iwvl], alpha[iwvl], 0.5)*dps[iwvl] #EE-50 en-squared diameter
            #ress[iwvl, 2]=moffat_ensquare_width(gamma[iwvl], alpha[iwvl], 0.8)*dps[iwvl] #EE-80 en-squared diameter
            ress[iwvl, 1]=moffat_slit_width(gamma[iwvl], alpha[iwvl], 0.5)*dps[iwvl] #EE-50 slit width
            ress[iwvl, 2]=moffat_slit_width(gamma[iwvl], alpha[iwvl], 0.8)*dps[iwvl] #EE-80 slit width
            #alpha is manually set to 3, so the ratio between fwhm and EE-p% is constant
            #not a good idea to fit over alpha
        np.savez(fn_cache, dps=dps, wvls=wvls, alpha=alpha, gamma=gamma, ress=ress)
    return ress,dps,wvls

def proc_psfs(fd, seeds=[1], **kargs):
    """
    Read PSFs from a folder and compute FWHM.
    Multiple seeds are averaged.
    """
    ressc=[]
    for seed in seeds:
        if not os.path.isfile(f"{fd}/Res_{seed}.done"):
            continue
        ress=[] #result array
        fr=[] #field point radius
        fx=[] #field point x
        fy=[] #field point y
        dps=None #PSF pixel sampling
        fns=natsorted(glob.glob(f"{fd}/evlpsfcl_{seed}_x*_y*.fits"))
        for fn in fns:
            m = re.search(r"_x([+-]?\d+(?:\.\d+)?)_y([+-]?\d+(?:\.\d+)?)\.fits$", fn)
            x, y = map(float, m.groups())
            fr.append(np.sqrt(x*x+y*y))
            fx.append(x)
            fy.append(y)
            res,dps,wvls=proc_psf(fn, **kargs)
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
        fn=f"{fd}/evlpsfol_{seed}.fits"
        datas=read(fn)
        ressol, dps, wvls=proc_psf(fn, **kargs)
    
        ressc.append({'avg':avg, 'std':err, 'cl':ress, 'ol':ressol,'fr':fr, 'fx':fx, 'fy':fy, 'dps':dps, 'wvls':wvls})
    if len(ressc)>1: #collect seeds
        ressc = {k: rms(np.array([d[k] for d in ressc if d is not None]), axis=0) for k in ressc[0]}
    elif len(ressc)==1:
        ressc = ressc[0]
    else:
        ressc = None
    return ressc
    
def strip_prefix(strings):
    if not strings:
        return strings

    prefix = os.path.commonprefix(strings)
    return [s[len(prefix):] for s in strings]
def longest_common_substring(strings):
    """Return the longest substring common to all strings."""
    if not strings:
        return ""

    shortest = min(strings, key=len)

    for length in range(len(shortest), 0, -1):
        for start in range(len(shortest) - length + 1):
            sub = shortest[start:start + length]
            if all(sub in s for s in strings):
                return sub
    return ""
def strip_common(strings):
    """
    Remove:
      1. common directory path
      2. common filename prefix
      3. common filename suffix

    Parameters
    ----------
    strings : list[str]

    Returns
    -------
    list[str]
    """
    if not strings:
        return []

    # Remove common directory path
    paths = [Path(s) for s in strings]
    parent = os.path.commonpath([str(p.parent) for p in paths])

    names = []
    for p in paths:
        try:
            rel = p.relative_to(parent)
        except ValueError:
            # Different drives (Windows), fall back to filename
            rel = p.name
        names.append(str(rel))

    # Remove common prefix
    prefix = os.path.commonprefix(names)
    names = [n[len(prefix):] for n in names]

    # Remove common suffix
    rev_prefix = os.path.commonprefix([n[::-1] for n in names])
    suffix = rev_prefix[::-1]

    if suffix:
        names = [n[:-len(suffix)] for n in names]

    # Remove longest common middle substring
    middle = longest_common_substring(names)
    if middle:
        names = [s.replace(middle, "...", 1) for s in names]

    return names    