"""
Measure the PSF size (FWHM and ensquared energy width) from MAOS PSF outputs.

MAOS saves the time averaged PSFs as FITS files: evlpsfcl_{seed}_x{x}_y{y}.fits
(closed loop, one file per evaluation direction at field location x, y in
arcsec) and evlpsfol_{seed}.fits (open loop, field averaged). Each file holds
one PSF per wavelength with header keywords WVL (wavelength in m), DP (PSF
pixel size in arcsec), and SUM (total PSF intensity).

Two methods are provided:

- proc_psf(): numerical. The FWHM is from the azimuthal average and the
  ensquared energy width is interpolated from the ensquared energy curve, both
  computed with FFT by libaos denc().
- proc_psf_fit(): fit a Moffat (or Gaussian + Moffat) model with moffat_fit()
  and compute the FWHM (moffat_fwhm) and ensquared energy width
  (moffat_square_width) from the model.

proc_psfs() processes all the directions of a simulation folder, averaging
PSFs of multiple seeds, and caches intermediate results next to the PSF files
(*.fits.npz) that are reused while newer than the PSF file.

Requires libaos (scripts/libaos.py with aolib.so from the MAOS build), numpy,
scipy and astropy. natsort is optional.
"""
import os
from pathlib import Path
import numpy as np
import glob


#from readbin import readbin as read
from libaos import read #this uses maos libaos.so
from astropy.modeling.models import Moffat2D
from astropy.modeling.fitting import TRFLSQFitter #handles bounded optim better than LevMarLSQFitter
from astropy.modeling import Fittable2DModel, Parameter
from scipy.special import betainc, erf
from scipy.integrate import dblquad
from scipy.optimize import brentq
from warnings import warn
import libaos as aos
try:
    from natsort import natsorted
except:
    natsorted=sorted
import re
def strip_prefix(strings):
    """
    Remove the common prefix, and then any leading underscores, from each string.

    Parameters
    ----------
    strings : list[str]

    Returns
    -------
    list[str]
        The input unchanged if it is empty.
    """
    if not strings:
        return strings

    prefix = os.path.commonprefix(strings)
    arr=[s[len(prefix):] for s in strings]
    arr = [s.lstrip("_") for s in arr]
    return arr
def strip_common(strings):
    """
    Shorten a list of file paths to the parts that differ, for use as labels.

    Removes, in order:
      1. the common directory path
      2. the common prefix of the remaining names
      3. the common suffix
      4. the longest substring common to all the names, replaced by "...",
         if it is longer than 5 characters

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
    if len(middle)>5:
        names = [s.replace(middle, "...", 1) for s in names]

    return names    
# import maos_client
def rms(*args, **kargs):
    """
    Root mean square of the input.

    Parameters
    ----------
    *args : array_like
        A single array, or several values that are combined into one array.
    **kargs
        Passed to np.mean, e.g. axis.

    Returns
    -------
    The rms, sqrt(mean(|x|^2)). Complex input is supported.
    """
    if len(args) > 1:
        arr = np.array(args)
    else:
        arr = np.array(args[0])
    return np.sqrt(np.mean(np.real(arr * np.conj(arr)), **kargs))
def fft_size(n, larger=True):
    """
    Return an even integer with no prime factors other than 2, 3, 5 and 7,
    which is efficient for FFT.

    Parameters
    ----------
    n : int
        The reference size.
    larger : bool
        True: the smallest such size >= n. False: the largest such size <= n
        (n must be at least 2).

    Returns
    -------
    int
    """
    def ok(m):
        for f in (2, 3, 5, 7):
            while m % f == 0:
                m //= f
        return m == 1
    m = max(int(np.ceil(n)) if larger else int(np.floor(n)), 2)
    m += m & 1 if larger else -(m & 1)
    while not ok(m):
        m += 2 if larger else -2
    return m
def auto_crop_roi(img, fraction=0.9, total=None):
    """
    Crop the image to the smallest square, centered on the FFT center, that
    encloses `fraction` of the total energy.

    The square of even size n has its pixel (n/2, n/2) at pixel (nx/2, ny/2)
    of the image (the FFT center, where MAOS puts the PSF). n is increased to
    a size convenient for FFT (fft_size()). If the energy is not reached, or n
    does not fit in the image, n is the largest such size that fits.

    Parameters
    ----------
    img : 2-D array
        The image, indexed as img[y, x].
    fraction : float
        Fraction of the total energy to enclose, in (0, 1].
    total : float, optional
        The total energy, e.g., the SUM keyword of a MAOS PSF, which includes
        the energy outside of the image. Defaults to the sum of img.

    Returns
    -------
    sub : 2-D array
        The cropped n*n image, a C-contiguous copy so that it can be modified
        in place, e.g. by libaos dshift2center(). If the image is smaller than
        2 pixels along either dimension, a copy of the whole image.
    (x0, y0) : tuple of int
        Position of sub[0, 0] in img.
    """
    ny, nx = img.shape
    if min(nx, ny) < 2:
        return np.array(img, order='C'), (0, 0)
    cy, cx = ny // 2, nx // 2
    #largest even half-size-symmetric square around (cy, cx) that fits
    nmax = 2 * min(cy, cx, ny - cy, nx - cx)
    if total is None:
        total = np.sum(img)
    #integral image: cs[i, j] is the sum of img[:i, :j]
    cs = np.zeros((ny + 1, nx + 1))
    cs[1:, 1:] = np.cumsum(np.cumsum(img, axis=0), axis=1)
    h = np.arange(1, nmax // 2 + 1) #half size
    ee = cs[cy + h, cx + h] - cs[cy - h, cx + h] - cs[cy + h, cx - h] + cs[cy - h, cx - h]
    ind = np.nonzero(ee >= fraction * total)[0]
    n = 2 * h[ind[0]] if ind.size else nmax
    n = fft_size(n)
    if n > nmax:
        n = fft_size(nmax, larger=False)
    y0, x0 = cy - n // 2, cx - n // 2

    return np.array(img[y0:y0+n, x0:x0+n], order='C'), (x0, y0)

class EllipticalMoffat2D(Fittable2DModel):
    """
    Elliptical Moffat with major axis and axis ratio parameterization:

        I(x, y) = amplitude * (1 + (x'/gamma)^2 + (y'/(q*gamma))^2)^(-alpha)

    where (x', y') are the coordinates relative to (x_0, y_0) rotated by theta
    so that x' is along the major axis.

    Parameters
    ----------
    amplitude : peak value
    x_0, y_0 : center (pixel)
    gamma : Moffat core width along the major axis (pixel)
    q : minor/major axis ratio, in (0, 1]
    alpha : Moffat power index. The total flux is finite only for alpha>1.
    theta : rotation of the major axis from the x axis (radian), in [-pi/2, pi/2]
    """
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
    """
    Sum of an elliptical Gaussian (core), an elliptical Moffat (halo) and a
    constant background. Both components share the center, axis ratio and
    rotation:

        I(x, y) = background
                + amplitude_g * exp(-0.5*((x'/sigma)^2 + (y'/(q*sigma))^2))
                + amplitude_m * (1 + (x'/gamma)^2 + (y'/(q*gamma))^2)^(-alpha)

    where (x', y') are the coordinates relative to (x_0, y_0) rotated by theta
    so that x' is along the major axis.

    Parameters
    ----------
    x_0, y_0 : center (pixel)
    amplitude_g, sigma : Gaussian peak and standard deviation along the major
        axis (pixel)
    amplitude_m, gamma, alpha : Moffat peak, core width along the major axis
        (pixel) and power index
    q : minor/major axis ratio
    theta : rotation of the major axis from the x axis (radian)
    background : constant background
    """
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
        x_0, y_0,
        amplitude_g, sigma,
        amplitude_m, gamma, alpha,
        q, theta, background
    ):#the parameters are passed in the order they are declared above
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

def moffat_fwhm(model):
    """
    Calculate the circularized FWHM of a PSF model.

    The half maximum radius along the major axis is found numerically from the
    radial profile (the peak is the sum of the component amplitudes). For
    elliptical models:

        FWHM = sqrt(FWHM_major * FWHM_minor) = FWHM_major * sqrt(q)

    Supported models (identified by their parameter names):
        Circular Moffat (astropy Moffat2D):
            amplitude, gamma, alpha
        Elliptical Moffat (EllipticalMoffat2D):
            amplitude, gamma, alpha, q
        Elliptical Gaussian + Moffat (EllipticalGaussianMoffat):
            amplitude_g, sigma, amplitude_m, gamma, alpha, q

    The background, center and rotation do not affect the result.

    Parameters
    ----------
    model : astropy model
        A fitted model, e.g., from moffat_fit().

    Returns
    -------
    float
        The FWHM in pixels.

    Raises
    ------
    TypeError
        For an unsupported model.
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
      
def moffat_slit_width(model, fraction=0.8):
    """
    Return the width of an infinitely long slit, centered on the PSF, that
    encloses `fraction` of the total PSF energy.

    The slit is aligned with the PSF major axis, so the width is along the
    minor axis. It is solved from the analytic slit enclosed energy of each
    component (a regularized incomplete beta function for the Moffat, erf for
    the Gaussian, weighted by their flux), in log(width) for robustness.

    Supported models are the same as for moffat_fwhm(). The background is not
    included.

    Parameters
    ----------
    model : astropy model
        A fitted model, e.g., from moffat_fit().
    fraction : float
        Enclosed energy fraction, in (0, 1).

    Returns
    -------
    float
        The slit width in pixels.

    Raises
    ------
    ValueError
        If fraction is outside (0, 1), alpha<=1 (infinite flux), or the total
        flux is not positive.
    TypeError
        For an unsupported model.
    RuntimeError
        If the width cannot be bracketed.
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
def moffat_flux(model):
    """
    Total flux of a fitted PSF model, integrated analytically to infinity.

    Moffat: pi*amplitude*gamma^2*q/(alpha-1). Gaussian: 2*pi*amplitude_g*sigma^2*q.
    Supported models are the same as for moffat_fwhm(). The background is not
    included.

    Parameters
    ----------
    model : astropy model
        A fitted model, e.g., from moffat_fit().

    Returns
    -------
    float
        The total flux in the units of the image times pixel^2.

    Raises
    ------
    ValueError
        If alpha<=1 (infinite flux).
    """
    names = model.param_names
    q = model.q.value if "q" in names else 1.0
    alpha = model.alpha.value
    if alpha <= 1:
        raise ValueError(f"alpha={alpha} <= 1. Moffat total flux is infinite.")
    if "amplitude_g" in names:
        return (2.0 * np.pi * model.amplitude_g.value * model.sigma.value**2 * q
            + np.pi * model.amplitude_m.value * model.gamma.value**2 * q / (alpha - 1.0))
    else:
        return np.pi * model.amplitude.value * model.gamma.value**2 * q / (alpha - 1.0)
def moffat_square_width(model, fraction=0.8):
    """
    Return the width of the square, aligned with the image axes and centered
    on the PSF, that encloses `fraction` of the total PSF energy.

    The enclosed energy is integrated numerically (dblquad) with the
    background set to 0 and normalized by moffat_flux(). The width is solved
    with brentq to a tolerance of 1e-4 times the larger of gamma and sigma.

    Supported models are the same as for moffat_flux().

    Parameters
    ----------
    model : astropy model
        A fitted model, e.g., from moffat_fit().
    fraction : float
        Enclosed energy fraction, in (0, 1).

    Returns
    -------
    float
        The square width in pixels.

    Raises
    ------
    ValueError
        If fraction is outside (0, 1) or alpha<=1.
    RuntimeError
        If the width cannot be bracketed.
    """
    if not 0 < fraction < 1:
        raise ValueError("fraction must be between 0 and 1")
    names = model.param_names
    params = [getattr(model, name).value for name in names]
    if "background" in names:
        params[names.index("background")] = 0
    x_0 = model.x_0.value
    y_0 = model.y_0.value
    total = moffat_flux(model)
    scale = model.gamma.value
    if "sigma" in names:
        scale = max(scale, model.sigma.value)

    def ee(width):
        half = width / 2
        num, _ = dblquad(
            lambda y, x: model.evaluate(x + x_0, y + y_0, *params),
            -half, half, -half, half,
        )
        return num / total

    w_hi = 2 * scale
    for _ in range(100):
        if ee(w_hi) >= fraction:
            break
        w_hi *= 2
    else:
        raise RuntimeError(f"Could not bracket fraction={fraction}.")
    return brentq(lambda w: ee(w) - fraction, 0, w_hi, xtol=1e-4 * scale)
def moffat_fit(img, alpha=None, name=None, shape=1):
    """
    Fit a Moffat model to a PSF image.

    The image is first cropped with auto_crop_roi() to speed up the fit. The
    model is fitted with bounded least squares (TRFLSQFitter), with alpha
    bounded to [1.5, 20]. A message is printed if the fit did not converge,
    the relative residual sum((img-model)^2)/sum(img^2) exceeds 0.05, or
    alpha approaches 1.

    Parameters
    ----------
    img : 2-D array
        The PSF image, indexed as img[y, x].
    alpha : float, optional
        Initial alpha for shape 0 and 1, and the fixed alpha for shape 2.
        Defaults to 3, 2.5 and 2 for shape 0, 1 and 2 respectively.
    name : str, optional
        Name, e.g. the file name, printed when the fit is not good.
    shape : int
        The model:
        0: circular Moffat (astropy Moffat2D),
        1: elliptical Moffat (EllipticalMoffat2D),
        2: elliptical Gaussian + Moffat (EllipticalGaussianMoffat) with alpha
           fixed.

    Returns
    -------
    astropy model
        The fitted model, with x_0 and y_0 in the pixel coordinates of img.
    """
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
    if alpha is None:
        alpha=(3, 2.5, 2)[shape]
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
        alpha=alpha,
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
            alpha=alpha,
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
    ierr= fitter.fit_info["status"]<1 #least_squares status: -1 or 0 is failure
    #Sanity check results
    model_img=fit(x,y)
    diff=np.sum((sub-model_img)**2)/np.sum(sub**2)
    if ierr or diff>0.05 or fit.alpha.value<=1+2e-7:
        print(f"Fitting is not good: {name}, diff={diff:.4g}, alpha={fit.alpha.value:.4g}, gamma={fit.gamma.value:.4g}")
        if ierr:
            print("nfev:", fitter.fit_info["nfev"])
            print("message:", fitter.fit_info["message"])
    
    # shift back to global coords
    fit.x_0.value += xoff
    fit.y_0.value += yoff
    
    return fit

def parse_header_float(headers, key):
    """
    Read a numerical keyword from a list of headers.

    Parameters
    ----------
    headers : list
        One header per PSF. Each is either a FITS header string made of
        80-character cards "KEY     = value / comment", or a dict mapping the
        key to a "value / comment" string.
    key : str
        The upper case keyword, e.g., 'DP', 'WVL' or 'SUM'.

    Returns
    -------
    np.ndarray
        The value from each header.

    Raises
    ------
    Exception
        If the value cannot be parsed (no '/' after the value).
    """
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
    """
    Read PSFs from one file, or average them over several files (e.g., seeds).

    Parameters
    ----------
    fn : str or list[str]
        The file name(s). All files must have the same wavelengths and PSF
        sizes.

    Returns
    -------
    datas : array
        One PSF per wavelength, as read by libaos read().
    headers : list
        The header of each PSF from the first file.
    """
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
def proc_psf_fit(fn, fn_cache=None, EEs=(0.5, 0.8), **kargs):
    """
    Read PSFs and compute their size by fitting a model with moffat_fit().

    The fitted models are cached in fn_cache and reused while the cache is
    newer than the (first) PSF file. The sizes are recomputed from the models
    every time.

    Parameters
    ----------
    fn : str or list[str]
        The PSF file, or a list of files to average (see average_psf()).
    fn_cache : str, optional
        The cache file (.npz). Defaults to fn+'.npz' for a single file. For a
        list of files no cache is used unless given.
    EEs : sequence of float
        Ensquared energy fractions.
    **kargs
        Passed to moffat_fit(), e.g., shape or alpha.

    Returns
    -------
    ress : ndarray, shape (nwvl, 1+len(EEs))
        Column 0 is the circularized FWHM (moffat_fwhm()). Column 1+i is the
        width of the square that encloses EEs[i] of the energy
        (moffat_square_width()). All in arcsec.
    dps : ndarray
        PSF pixel size (arcsec) for each wavelength.
    sums : ndarray
        Total PSF intensity for each wavelength.
    wvls : ndarray
        Wavelengths (micron).
    model : ndarray of object
        The fitted model for each wavelength.
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
            model[iwvl]=moffat_fit(datas[iwvl], name=fn0, **kargs)
        if fn_cache is not None:
            np.savez(fn_cache, dps=dps, sums=sums, wvls=wvls, model=model)
    #no longer caching ress
    nwvl=wvls.size
    ress=np.zeros((nwvl,1+len(EEs)))
    for iwvl in range(nwvl):
        ress[iwvl, 0]=moffat_fwhm(model[iwvl])*dps[iwvl] #FWHM
        #ress[iwvl, 1]=moffat_encircle_width(gamma[iwvl], alpha[iwvl], 0.5)*dps[iwvl] #EE-50 diameter
        #ress[iwvl, 2]=moffat_encircle_width(gamma[iwvl], alpha[iwvl], 0.8)*dps[iwvl] #EE-80 en-circled diameter
        for iEE in range(len(EEs)):
            ress[iwvl, 1+iEE]=moffat_square_width(model[iwvl], EEs[iEE])*dps[iwvl] #EE- en-squared width
            #ress[iwvl, 2+iEE*2]=moffat_slit_width(model[iwvl], EEs[iEE])*dps[iwvl] #EE- slit width
        
    return ress,dps,sums,wvls,model
def proc_psf(fn, fn_cache=None, EEs=(0.5, 0.8), **kargs):
    """
    Read PSFs and compute their size numerically with libaos denc().

    For each wavelength the PSF is shifted to put its center of gravity at the
    FFT center, then for r=0, 1, ..., nx/2-1 pixels:
      - azavg: the azimuthal average at radius r, normalized by the center value;
      - ensquare: the energy within a square of width 2r, normalized by the
        total PSF intensity (header SUM), so the energy scattered outside of
        the PSF image is accounted for.
    These are cached in fn_cache and reused while the cache is newer than the
    (first) PSF file. The sizes are interpolated from them every time.

    Parameters
    ----------
    fn : str or list[str]
        The PSF file, or a list of files to average (see average_psf()).
    fn_cache : str, optional
        The cache file (.npz). Defaults to fn+'.npz' for a single file. For a
        list of files no cache is used unless given.
    EEs : sequence of float
        Ensquared energy fractions.
    **kargs
        Not used. Accepted so that proc_psfs() can pass the same arguments to
        proc_psf() and proc_psf_fit().

    Returns
    -------
    ress : ndarray, shape (nwvl, 1+len(EEs))
        Column 0 is the FWHM, twice the radius where azavg drops to 0.5.
        Column 1+i is the width of the square that encloses EEs[i] of the
        energy. All in arcsec.
    dps : ndarray
        PSF pixel size (arcsec) for each wavelength.
    sums : ndarray
        Total PSF intensity for each wavelength.
    wvls : ndarray
        Wavelengths (micron).
    enc : ndarray of object
        For each wavelength a dict with 'r' (pixel), 'azavg', 'ensquare' and
        'enslit' (currently None).
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
        datas, headers=average_psf(fn)
        dps=parse_header_float(headers, 'DP') #pixel width
        wvls=parse_header_float(headers, 'WVL')*1e6 #wavelength, convert to micron
        sums=parse_header_float(headers, 'SUM') #total intensity
        nwvl=datas.shape[0]
        enc=np.empty((nwvl), dtype=object)
        
        for iwvl in range(nwvl):
            img, _=auto_crop_roi(datas[iwvl], total=sums[iwvl]) #enclose 90% of the total energy
            nx, ny=img.shape
            if nx!=datas[iwvl].shape[0]:
                print(f'size reduced from {datas[iwvl].shape[0]} to {nx}')
            aos.dshift2center(img, 0.5, 0.5) #make CoG at FFT zero frequency
            r=np.arange(0, nx/2, 1)
            azavg=aos.denc(img, r, -1, 0)
            azavg/=azavg[0]
            
            #enslit=(aos.denc(img, r*2, 2, 0)+aos.denc(img.T, r*2, 2, 0))/(2*sums[iwvl])
            enslit=None
            ensquare=aos.denc(img, r*2, 0, 0)/sums[iwvl];
            enc[iwvl]={'r':r, 'azavg':azavg, 'ensquare':ensquare, 'enslit':enslit}
        if fn_cache is not None:
            np.savez(fn_cache, dps=dps, sums=sums, wvls=wvls, enc=enc)
    #no longer caching ress
    nwvl=wvls.size
    ress=np.zeros((nwvl,1+len(EEs)))
    for iwvl in range(nwvl):
        azavg=enc[iwvl]['azavg']
        r=enc[iwvl]['r']
        ensquare=enc[iwvl]['ensquare']
        #enslit=enc[iwvl]['enslit']
        #np.interp require data to be ascending
        ress[iwvl, 0]=np.interp(0.5, azavg[::-1], r[::-1])*2*dps[iwvl] #FWHM
        for iEE in range(len(EEs)):
            ress[iwvl, 1+iEE]=np.interp(EEs[iEE], ensquare, r)*2*dps[iwvl] #Ensquared width
            #ress[iwvl, 2+iEE*2]=np.interp(EEs[iEE], enslit, r)*2*dps[iwvl] #Ensqlited width

    return ress,dps,sums, wvls,enc

def proc_psfs(fd, seeds=[1], use_fit=0, fdol=None,**kargs):
    """
    Compute the PSF size of all the evaluation directions in a MAOS output
    folder.

    The closed loop PSFs evlpsfcl_{seed}_x{x}_y{y}.fits are found in fd for the
    first seed. For each direction, and for the open loop PSF
    evlpsfol_{seed}.fits in fdol, the PSFs of all seeds are averaged before
    being processed by proc_psf() (use_fit=0) or proc_psf_fit() (use_fit=1).
    The intermediate results are cached in
    {fd}/evlpsfcl_{seeds}_x{x}_y{y}_{enc|moffat}.fits.npz and
    {fdol}/evlpsfol_{seeds}_{enc|moffat}.fits.npz, where {seeds} are the
    seeds joined by '_'.

    Parameters
    ----------
    fd : str
        The folder with the closed loop PSFs.
    seeds : list[int]
        The seeds to average.
    use_fit : int
        0: numerical method (proc_psf()). 1: model fitting (proc_psf_fit()).
    fdol : str, optional
        The folder with the open loop PSF. Defaults to fd.
    **kargs
        Passed to proc_psf() or proc_psf_fit(), e.g., EEs.

    Returns
    -------
    dict or None
        None if no closed loop PSF is found. Otherwise, with nfield
        directions sorted by x then y, and nres=1+len(EEs) columns as returned
        by proc_psf():
        'cl' : (nfield, nwvl, nres) closed loop results (arcsec)
        'avg', 'std' : (nwvl, nres) mean and standard deviation of 'cl'
            over the field
        'ol' : (nwvl, nres) open loop results (arcsec)
        'fx', 'fy', 'fr' : (nfield,) field location x, y and radius (arcsec)
        'dps' : PSF pixel size (arcsec) for each wavelength
        'wvls' : wavelengths (micron)
    """
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
    return {'avg':avg, 'std':err, 'cl':ress, 'ol':ressol,'fr':fr, 'fx':fx, 'fy':fy, 'dps':dps, 'wvls':wvls}
    
def concat_dict(ress):
    """
    Stack a list of dicts with the same keys (e.g., from proc_psfs() for
    several folders) into a dict of arrays with a new leading axis. None
    entries are skipped.
    """
    return {k: np.array([d[k] for d in ress if d is not None]) for k in ress[0]}
