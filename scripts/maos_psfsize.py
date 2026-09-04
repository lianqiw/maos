import os
from pathlib import Path
import numpy as np
import glob


#from readbin import readbin as read
from libaos import read #this uses maos libaos.so
from astropy.modeling.models import Moffat2D, Const2D
from astropy.modeling.fitting import LevMarLSQFitter #handles bounded optim poorly
from astropy.modeling.fitting import TRFLSQFitter #handles bounded optim better
from astropy.modeling import Fittable2DModel, Parameter
from scipy.special import betaincinv, gamma as gammafun, betainc
from scipy.integrate import dblquad, quad
from scipy.optimize import brentq
from cycler import cycler
from warnings import warn
import libaos as aos
try:
    from natsort import natsorted
except:
    natsorted=sorted
import re
def strip_prefix(strings):
    if not strings:
        return strings

    prefix = os.path.commonprefix(strings)
    arr=[s[len(prefix):] for s in strings]
    arr = [s.lstrip("_") for s in arr]
    return arr
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
# import maos_client
def rms(*args, **kargs):
    """compute rms of input"""
    if len(args) > 1:
        arr = np.array(args)
    else:
        arr = np.array(args[0])
    return np.sqrt(np.mean(np.real(arr * np.conj(arr)), **kargs))
def auto_crop_roi(img, threshold_rel=0.1, max_radius=None):
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
class EllipticalGaussianMoffat(Fittable2DModel):
    x_0 = Parameter()
    y_0 = Parameter()

    amplitude_g = Parameter()
    sigma = Parameter()

    amplitude_m = Parameter()
    gamma = Parameter()
    alpha = Parameter()

    q = Parameter()

    theta = Parameter()
    background = Parameter()

    @staticmethod
    def evaluate(
        x, y,
        amplitude_g, sigma,
        amplitude_m, gamma, alpha,
        q, x_0, y_0, theta, background
    ):
        dx = x - x_0
        dy = y - y_0

        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)

        # Major/minor-axis coordinates
        xp =  dx * cos_theta + dy * sin_theta
        yp = -dx * sin_theta + dy * cos_theta

        # Elliptical Gaussian
        gaussian = amplitude_g * np.exp(
            -0.5 * (
                (xp / sigma)**2 +
                (yp / (q * sigma))**2
            )
        )

        # Elliptical Moffat
        moffat = amplitude_m * (
            1
            + (xp / gamma)**2
            + (yp / (q * gamma))**2
        )**(-alpha)

        return background + gaussian + moffat

def global_fwhm(model):
    """
    Calculate the circularized FWHM of a PSF model.

    Supported models
    ----------------
    Circular Moffat:
        amplitude, gamma, alpha

    Elliptical Moffat:
        amplitude, gamma, alpha, q

    Gaussian + Moffat:
        amplitude_g, sigma,
        amplitude_m, gamma, alpha, q

    For elliptical models:
        FWHM_global = sqrt(FWHM_major * FWHM_minor)
                     = FWHM_major * sqrt(q)

    Assumes the PSF is centered and the elliptical axes are
    represented by the q parameter.
    """

    names = model.param_names

    # ---------------------------------------------------------
    # Gaussian + Moffat
    # ---------------------------------------------------------
    if all(name in names for name in (
        "amplitude_g",
        "sigma",
        "amplitude_m",
        "gamma",
        "alpha",
    )):

        A_g = model.amplitude_g.value
        sigma = model.sigma.value

        A_m = model.amplitude_m.value
        gamma = model.gamma.value
        alpha = model.alpha.value

        q = model.q.value if "q" in names else 1.0

        peak = A_g + A_m

        def profile(r):
            gaussian = A_g * np.exp(
                -0.5 * (r / sigma)**2
            )

            moffat = A_m * (
                1 + (r / gamma)**2
            )**(-alpha)

            return gaussian + moffat

    # ---------------------------------------------------------
    # Moffat
    # ---------------------------------------------------------
    elif all(name in names for name in (
        "amplitude",
        "gamma",
        "alpha",
    )):

        A = model.amplitude.value
        gamma = model.gamma.value
        alpha = model.alpha.value

        q = model.q.value if "q" in names else 1.0

        peak = A

        def profile(r):
            return A * (
                1 + (r / gamma)**2
            )**(-alpha)

    else:
        raise TypeError(
            "Unsupported model. Expected a circular/elliptical "
            "Moffat or Gaussian + Moffat model."
        )

    # ---------------------------------------------------------
    # Find half-maximum radius along major axis
    # ---------------------------------------------------------

    half = peak / 2

    # Automatically find an upper bracket
    r_hi = 1.0

    while profile(r_hi) > half:
        r_hi *= 2

        if r_hi > 1e10:
            raise RuntimeError(
                "Could not bracket the half-maximum."
            )

    r_half = brentq(
        lambda r: profile(r) - half,
        0,
        r_hi,
    )

    fwhm_major = 2 * r_half

    # Circularized FWHM
    return fwhm_major * np.sqrt(q)
      
def slit_width(model, fraction=0.8):
    """
    Return the width of an infinitely long slit that encloses
    `fraction` of the total PSF energy.

    The slit is aligned with the PSF major axis.

    Supported models:
        1. Circular Moffat:
           amplitude, gamma, alpha

        2. Elliptical Moffat:
           amplitude, gamma, alpha, q

        3. Elliptical Gaussian + Moffat:
           amplitude_g, sigma,
           amplitude_m, gamma, alpha, q
    """

    if not 0 < fraction < 1:
        raise ValueError("fraction must be between 0 and 1")

    names = model.param_names

    # ---------------------------------------------------------
    # Moffat EEF for an infinitely long slit
    # ---------------------------------------------------------

    def moffat_eef(log_t, alpha):
        """
        t = W / (2 q gamma)

        EEF = 1 - I_s(alpha-1, 1/2)

        where
            s = 1 / (1 + t^2)

        This form is numerically more stable as alpha -> 1.
        """
        log_s = -np.logaddexp(0.0, 2.0 * log_t)
        s = np.exp(log_s)

        return 1.0 - betainc(
            alpha - 1.0,
            0.5,
            s,
        )

    # ---------------------------------------------------------
    # Circular / elliptical Moffat
    # ---------------------------------------------------------

    if all(p in names for p in (
        "amplitude",
        "gamma",
        "alpha",
    )):

        gamma_m = model.gamma.value
        alpha = model.alpha.value

        q = model.q.value if "q" in names else 1.0

        if alpha <= 1:
            raise ValueError(
                f"alpha={alpha} <= 1. "
                "Moffat total flux is infinite."
            )

        def eef_log_width(log_w):

            log_t = (
                log_w
                - np.log(2.0 * q * gamma_m)
            )

            return moffat_eef(
                log_t,
                alpha,
            )

        scale = 2.0 * q * gamma_m

    # ---------------------------------------------------------
    # Gaussian + Moffat
    # ---------------------------------------------------------

    elif all(p in names for p in (
        "amplitude_g",
        "sigma",
        "amplitude_m",
        "gamma",
        "alpha",
    )):

        A_g = model.amplitude_g.value
        sigma = model.sigma.value

        A_m = model.amplitude_m.value
        gamma_m = model.gamma.value
        alpha = model.alpha.value

        q = model.q.value if "q" in names else 1.0

        if alpha <= 1:
            raise ValueError(
                f"alpha={alpha} <= 1. "
                "Moffat total flux is infinite."
            )

        # Integrated 2-D flux of each component
        flux_g = (
            2.0 * np.pi
            * A_g
            * sigma**2
            * q
        )

        flux_m = (
            np.pi
            * A_m
            * gamma_m**2
            * q
            / (alpha - 1.0)
        )

        flux_total = flux_g + flux_m

        if flux_total <= 0:
            raise ValueError(
                "Total PSF flux must be positive."
            )

        weight_g = flux_g / flux_total
        weight_m = flux_m / flux_total

        def eef_log_width(log_w):

            width = np.exp(log_w)

            # Gaussian slit EEF
            eef_g = erf(
                width
                / (
                    2.0
                    * np.sqrt(2.0)
                    * q
                    * sigma
                )
            )

            # Moffat slit EEF
            log_t = (
                log_w
                - np.log(2.0 * q * gamma_m)
            )

            eef_m = moffat_eef(
                log_t,
                alpha,
            )

            return (
                weight_g * eef_g
                + weight_m * eef_m
            )

        scale = 2.0 * q * max(
            sigma,
            gamma_m,
        )

    else:
        raise TypeError(
            "Unsupported model."
        )

    # ---------------------------------------------------------
    # Solve EEF(W) = fraction in log(W)
    # ---------------------------------------------------------

    log_w_lo = np.log(scale) - 10.0
    log_w_hi = np.log(scale)

    for _ in range(1000):

        if eef_log_width(log_w_hi) >= fraction:
            break

        log_w_hi += np.log(2.0)

    else:
        raise RuntimeError(
            f"Could not bracket fraction={fraction}."
        )

    log_w = brentq(
        lambda log_w:
            eef_log_width(log_w) - fraction,
        log_w_lo,
        log_w_hi,
    )

    return np.exp(log_w)
def fit_moffat_roi(img, alpha=3, name=None, shape=1):
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
    if shape==0: #axial symmetric
        model = Moffat2D(
            amplitude=sub.max() ,#- B0,
            x_0=ix,
            y_0=iy,
            gamma=min(nx, ny) / 4,
            alpha=alpha
        ) #+ Const2D(amplitude=B0)
    elif shape==1: #elliptical
        model = EllipticalMoffat2D(
        amplitude=np.max(sub),
        x_0=ix,
        y_0=iy,
        gamma=3.0,
        q=0.8,
        alpha=2.5,
        theta=0.0,
    )
    elif shape==2: #Elliptical moffat + ellitical Gaussian
        model = EllipticalGaussianMoffat(
            x_0=ix,
            y_0=iy,
            #for Gaussian
            amplitude_g=np.max(sub)/2,
            sigma=2,
            #for Moffat
            amplitude_m=np.max(sub)/2,
            gamma=5,
            alpha=2,
            #Shared ellipticity. q is minor/major ratio
            q=0.8,
            theta=0,
            background=0,
        )
        model.alpha.fixed = True
    # bounded optimization
    model.alpha.bounds = (1.5, 20) #below 1 not integratable, >>10 degenerate with gamma
    fitter = TRFLSQFitter()
    #weights=(sub/np.max(sub))
    #weights[weights<.01]=0
    fit = fitter(model, x, y, sub)#, weights=weights)
    ierr= fitter.fit_info["status"]<1
    #Sanity check results
    model_img=fit(x,y)
    diff=np.sum((sub-model_img)**2)/np.sum(sub**2)
    if ierr or diff>0.05 or fit.alpha.value<=1+2e-7:
        print(f"Fitting is not good: {name}, diff={diff:.4g}, alpha={fit.alpha.value:.4g}, gamma={fit.gamma.value:.4g}")
        if ierr>4:
            print("nfev:", fitter.fit_info["nfev"])
            print("message:", fitter.fit_info["message"])
    
    # shift back to global coords
    fit.x_0.value += xoff
    fit.y_0.value += yoff
    
    return fit

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
            warn(f'dp: {dp}')
            raise(Exception('Unable to parse header'))
        dp=float(dp[:dp.find('/')])
        res.append(dp)
    return np.stack(res)
def average_psf(fn):
    if isinstance(fn, str):
        datas=read(fn)
    else:
        datas=read(fn[0])
    headers=[data.header for data in datas]
    if not isinstance(fn, str):
        count=1
        for fn2 in fn[1:]:
            datas+=read(fn2)
            count+=1
        if count>1:
            datas/=count

    return datas, headers
def proc_psf_fit(fn, fn_cache=None, **kargs):
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
            data=np.load(fn_cache, allow_pickle=True)
            dps=data["dps"]
            wvls=data["wvls"]
            sums=data["sums"]
            model=data["model"]
            skip=1
        except Exception as e:
            print(f'Loading cached results from {fn_cache} failed: {e}')
            pass

    if skip==0:
        datas, headers=average_psf(fn)
        dps=parse_header_float(headers, 'DP') #pixel width
        wvls=parse_header_float(headers, 'WVL')*1e6 #wavelength, convert to micron
        sums=parse_header_float(headers, 'SUM') #total intensity
        nwvl=datas.shape[0]
        model=np.empty((nwvl), dtype=object)
        
        for iwvl in range(nwvl):
            #res=calc_fwhm_gaussian(datas[iwvl], dx=dps[iwvl])
            #res=maos_utils.print_psf_metrics(directory=f"{fd}/", x=-90, y=0, ee=50, seed=1)
            model[iwvl]=fit_moffat_roi(datas[iwvl], name=fn0, **kargs)
        if fn_cache is not None:
            np.savez(fn_cache, dps=dps, sums=sums, wvls=wvls, model=model)
    #no longer caching ress
    nwvl=wvls.size
    ress=np.zeros((nwvl,3))
    for iwvl in range(nwvl):
        ress[iwvl, 0]=global_fwhm(model[iwvl])*dps[iwvl] #FWHM
        #ress[iwvl, 1]=moffat_encircle_width(gamma[iwvl], alpha[iwvl], 0.5)*dps[iwvl] #EE-50 diameter
        #ress[iwvl, 2]=moffat_encircle_width(gamma[iwvl], alpha[iwvl], 0.8)*dps[iwvl] #EE-80 en-circled diameter
        #ress[iwvl, 1]=moffat_ensquare_width(gamma[iwvl], alpha[iwvl], 0.5)*dps[iwvl] #EE-50 en-squared diameter
        #ress[iwvl, 2]=moffat_ensquare_width(gamma[iwvl], alpha[iwvl], 0.8)*dps[iwvl] #EE-80 en-squared diameter
        ress[iwvl, 1]=slit_width(model[iwvl], 0.5)*dps[iwvl] #EE-50 slit width
        ress[iwvl, 2]=slit_width(model[iwvl], 0.8)*dps[iwvl] #EE-80 slit width
        #alpha is manually set to 3, so the ratio between fwhm and EE-p% is constant
        #not a good idea to fit over alpha
    
    #print(f'Save results to {fn_cache}')
    return ress,dps,sums,wvls,model
def proc_psf(fn, fn_cache=None, **kargs):
    """
        Read PSF from a file and compute FWHM and ENC using FFT
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
            data=np.load(fn_cache, allow_pickle=True)
            dps=data["dps"] #pixel scale
            sums=data["sums"] #total intensity
            wvls=data["wvls"] #wavelength 
            enc=data["enc"] #processed data
            r=enc[0]['r'] #check
            skip=1
        except Exception as e:
            print(f'Loading cached results from {fn_cache} failed: {e}')
            pass

    if skip==0:
        try:
            datas, headers=average_psf(fn)
            dps=parse_header_float(headers, 'DP') #pixel width
            wvls=parse_header_float(headers, 'WVL')*1e6 #wavelength, convert to micron
            sums=parse_header_float(headers, 'SUM') #total intensity
            nwvl=datas.shape[0]
            enc=np.empty((nwvl), dtype=object)
            
            for iwvl in range(nwvl):
                img=datas[iwvl]
                aos.dshift2center(img, 0.5, 0.5) #make CoG at FFT zero frequency
                nx, ny=img.shape
                r=np.arange(0, nx/2, 0.5)
                azavg=aos.denc(img, r, -1, 0)
                azavg/=azavg[0]
                
                enslit=(aos.denc(img, r*2, 2, 0)+aos.denc(img.T, r*2, 2, 0))/(2*sums[iwvl])
                ensquare=aos.denc(img, r*2, 0, 0)/sums[iwvl];
                enc[iwvl]={'r':r, 'azavg':azavg, 'ensquare':ensquare, 'enslit':enslit}
            if fn_cache is not None:
                np.savez(fn_cache, dps=dps, sums=sums, wvls=wvls, enc=enc)
        except Exception as e:
            print(f'Process {fn} failed: {e}')
            return None
    #no longer caching ress
    nwvl=wvls.size
    ress=np.zeros((nwvl,3))
    for iwvl in range(nwvl):
        azavg=enc[iwvl]['azavg']
        r=enc[iwvl]['r']
        ensquare=enc[iwvl]['ensquare']
        enslit=enc[iwvl]['enslit']
        #np.interp require data to be ascending
        ress[iwvl, 0]=np.interp(0.5, azavg[::-1], r[::-1])*2*dps[iwvl] #FWHM
        ress[iwvl, 1]=np.interp(0.8, ensquare, r)*2*dps[iwvl] #Ensquared 80% width
        ress[iwvl, 2]=np.interp(0.8, enslit, r)*2*dps[iwvl] #Ensqlited 80% width
    return ress,dps,sums, wvls,enc

def proc_psfs(fd, seeds=[1], use_fit=0, fdol=None,**kargs):
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
        fn_cache=f"{fd}/evlpsfcl_{'_'.join(map(str, seeds))}_x{x}_y{y}{'_moffat' if use_fit else '_enc'}.fits.npz"
        if use_fit:
            tmp=proc_psf_fit(fns2, fn_cache, **kargs)
        else:
            tmp=proc_psf(fns2, fn_cache, **kargs)
        if tmp is not None:
            res,dps,_,wvls,_=tmp
            ress.append(res)
    if len(ress)==0:
        return None
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
    if fdol is None:
        fdol=fd
    fn=f"{fdol}/evlpsfol_{seeds[0]}.fits"
    fns2=[fn]
    for seed in seeds[1:]:
        fns2.append(f"{fdol}/evlpsfol_{seed}.fits")
    fn_cache=f"{fdol}/evlpsfol_{'_'.join(map(str, seeds))}{'_moffat' if use_fit else '_enc'}.fits.npz"
    if use_fit:
        ressol, dps, _, wvls, _=proc_psf_fit(fns2, fn_cache, **kargs)
    else:
        ressol, dps, _, wvls, _=proc_psf(fns2, fn_cache, **kargs)
    ressc.append({'avg':avg, 'std':err, 'cl':ress, 'ol':ressol,'fr':fr, 'fx':fx, 'fy':fy, 'dps':dps, 'wvls':wvls})
    if len(ressc)>1: #collect seeds
        ressc = {k: rms(np.array([d[k] for d in ressc if d is not None]), axis=0) for k in ressc[0]}
    elif len(ressc)==1:
        ressc = ressc[0]
    else:
        ressc = None
    return ressc
    
