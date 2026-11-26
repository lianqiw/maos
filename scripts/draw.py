#!/usr/bin/env python3
"""
This module contains routines to embed opds defined on coordinates to 2-d array and plot.

Usage:
    from draw import draw
    draw(img)              #plot a 2-d array (a 1-d array of n*n is reshaped to n*n)
    draw(loc)              #plot the points of loc (2*n array of x and y coordinates)
    draw(loc, opd)         #embed opd (defined on loc) on a 2-d grid and plot
    draw(cell)             #plot each element of a cell (object) array in subplots (up to 60)
    python3 draw.py file1.bin [file2.bin ...]  #plot files read with readbin

Keyword arguments of draw():
    ext=[x0,x1,y0,y1]: extent of the image (passed to imshow).
    dx=d: sampling of the image (or of loc) to compute the extent.
    center=n or (ny,nx): crop (or pad) the image to this size around the center.
    clim=(min,max): color limits.
    colorbar=0: disable the colorbar.
    nrow=, ncol= (or nx=): layout of the subplots for multiple images.
    figsize=: figure size for new figures.

Requirements: numpy, scipy, matplotlib.
"""
#from aolib import *
import sys
import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse
#from aotools import center #cyclic import
#import struct

def iscell(arr):
    """Whether arr is a cell array (numpy object array)."""
    return isinstance(arr, np.ndarray) and arr.dtype.name=='object'
    
def coord2grid(x, **kargs):
    """Map coordinates x to integer indices on a regular grid.

    Args:
        x: 1-d array of coordinates.
        dx (keyword): grid sampling. If absent or 0, it is estimated as the
            median of the spacing between sorted unique coordinates.

    Returns:
        (ix, nx, dx, xi): integer index of each coordinate, number of grid
        points, grid sampling, and [min, max] of x.
    """
    xmin = x.min()
    xmax = x.max()
    x2 = x-xmin
    if x2.max() > 0:
        if 'dx' in kargs and kargs['dx']!=0:
            dx=kargs['dx']
        else:
            tol=(xmax-xmin)*1e-6
            di=np.diff(np.sort(x2))
            dx=np.median(di[di>tol])
            #dx2=x2[x2>0].min() #assume regular
            #dx=(xmax-xmin)*np.sqrt(np.pi*0.5/x.size) #assume evenly distributed
            #if(abs(dx-dx2)<(dx+dx2)*0.1):
            #    dx=dx2
        dxi = 1./dx
        nx = np.round((xmax-xmin)*dxi+1.).astype(int)
        ix = np.floor(x2*dxi+1e-5).astype(int)
    else:
        ix = np.zeros(x.shape, dtype=int)
        nx = 1
        dx = 1
    #xi = [round(xmin/dx)*dx, round(xmax/dx)*dx]
    xi = [xmin, xmax]
    return (ix, nx, dx, xi)


def unwrap(arg):
    """Squeeze singleton dimensions and unwrap single element cells."""
    while isinstance(arg, np.ndarray):
        arg=np.squeeze(arg)
        if iscell(arg) and arg.size==1:
            arg=arg.flat[0]
        else:
            break
    return arg

def isvec(arg):
    """Whether arg is a vector (1-d, or 2-d with a singleton dimension)."""
    return arg.ndim == 1 or (arg.ndim == 2 and (arg.shape[0] == 1 or arg.shape[1] == 1))

def center(img, n):
    """Crop (or zero pad) img to n*n (or n[0]*n[1]) pixels around its center."""
    if np.isscalar(n):
        n=(n, n)
    out=np.zeros(n, dtype=img.dtype)
    ny, nx=img.shape
    oy, ox=(ny-n[0])//2, (nx-n[1])//2 #offset of out in img
    y0, x0=max(oy, 0), max(ox, 0)
    y1, x1=min(oy+n[0], ny), min(ox+n[1], nx)
    out[y0-oy:y1-oy, x0-ox:x1-ox]=img[y0:y1, x0:x1]
    return out


def isloc(arg):
    """Whether arg looks like a loc: a 2*n (or n*2) numeric array (n>2) of
    coordinates that are mostly distinct points on a regular grid."""
    if not isinstance(arg, np.ndarray):
        return False
    if iscell(arg):
        return False
    if arg.ndim !=2:
        return False
    if arg.shape[0] == 2 and arg.shape[1]>2:
        loc=arg
    elif arg.shape[0]>2 and arg.shape[1]==2:
        loc=arg.T
    else:
        return False
    (ix, nx, dx, xi) = coord2grid(loc[0])
    (iy, ny, dy, yi) = coord2grid(loc[1])
    
    nelem=np.unique(ix+iy*nx).size
    if nelem*10<loc.shape[1]:
        return False
    
    return True

def locembed(loc, opd, return_ext=0, **kargs):
    """Embed opd defined on coordinates loc onto a 2-d grid (used by draw(loc, opd)).

    Args:
        loc: 2*n (or n*2) array of x and y coordinates.
        opd: values on loc: a vector of n, an array of n*m or m*n (m frames),
            or a cell array of them. Sparse arrays and lists are accepted.
        return_ext: also return the extent.
        dx (keyword): grid sampling (default: estimated from loc).

    Returns:
        ims (2-d array of ny*nx with NaN where loc has no point, or an array
        of them for several frames or cells), and ext=[x0, x1, y0, y1] (edges
        of the pixels) if return_ext is set. None if the size of opd is not a
        multiple of n.
    """
    # draw(loc, opd): embed opd onto loc
    if loc.shape[0]>2 and loc.shape[1]==2: #need a transpose
        loc=loc.T
    if scipy.sparse.issparse(opd):
        opd=opd.toarray()
    elif isinstance(opd, list):
        opd=np.array(opd)
    nloc = loc.shape[1]
    if iscell(opd): #cell
        ims = []
        ext = None
        for opdi in opd:
            if opdi is not None and np.size(opdi)>0:
                ims_i,ext=locembed(loc, opdi, 1, **kargs)
                ims.append(ims_i)
        ims=np.asarray(ims)
        if ims.size == 1:
            ims=ims[0]
        if return_ext:
            return ims, ext
        else:
            return ims
    elif opd.ndim==2 and (opd.shape[0]==nloc or opd.shape[0]==nloc*2) and opd.shape[1]!=nloc:
        opd=opd.T
    if not opd.size % nloc == 0:
        print("data has wrong shape", opd.shape)
        return None
    if not opd.flags['FORC']:
        opd=np.ascontiguousarray(opd)
    nframe = int(opd.size/nloc)
    opd=opd.reshape(nframe, nloc) #reshape() may do copy. assign .shape prevents copy

    # print('opd=',opd.shape)
    ims = np.empty(nframe, dtype=object)
    (ix, nx, dx, xi) = coord2grid(loc[0], **kargs)
    (iy, ny, dy, yi) = coord2grid(loc[1], **kargs)
    
    for iframe in range(nframe):
        im = np.full((nx*ny),np.nan)
        im[ix+iy*nx] = opd[iframe, :]
        im.shape = (ny, nx)  # fixed 2020-10-09
        ims[iframe] = im
    ext = np.r_[xi[0]-dx/2, xi[1]+dx/2, yi[0]-dy/2, yi[1]+dy/2]
    if nframe == 1:
        ims=ims[0]
    if return_ext:
        return ims, ext
    else:
        return ims
def isimg(arg0):
    """Whether arg0 can be drawn as a single image (or curve): a sparse array, a
    1-d numeric array, or a 2-d numeric array whose aspect ratio is within 100."""
    if scipy.sparse.issparse(arg0):
        return True
    elif isinstance(arg0, np.ndarray):
        while (iscell(arg0) or arg0.ndim>2) and arg0.shape[0]==1:
            arg0=arg0[0]
        return iscell(arg0)==False and (arg0.ndim==1 or (arg0.ndim==2 and arg0.shape[0]*100>arg0.shape[1] and arg0.shape[1]*100>arg0.shape[0]))
def draw(*args, **kargs):
    """Plot arrays, cell arrays, or OPDs defined on coordinates (see the module docstring).

    Returns:
        draw(loc, opd): (ims, ext) as returned by locembed.
        draw(img): the image (or curve) drawn.
        draw(cell): the cell array drawn.
        None if nothing is drawn.
    """
    if(len(args)==0):
        return
    #if not 'keep' in kargs or kargs['keep'] == 0:
    #    plt.clf()
    arg0=unwrap(args[0])
    if isloc(arg0):  # first argument is loc
        loc=arg0
        if loc.shape[0]>2 and loc.shape[1]==2: #need a transpose
            loc=loc.T
        if len(args) == 1:
            # plot grid
            plt.plot(loc[0], loc[1], '+')
            #plt.axis('scaled')  # better than equal
            plt.gca().axis('tight')
            plt.gca().set_aspect('equal', adjustable='box')
            plt.xlabel('x (m)')
            plt.ylabel('y (m)')
            return None
        elif len(args) >= 2:
            ct=0
            ims, ext = None, None
            for arg1 in args[1:]:
                ct+=1
                if ct>1:
                    plt.figure(figsize=kargs.get('figsize'))
                res = locembed(loc, arg1, return_ext=1, **kargs)
                if res is None: #wrong shape. message is printed by locembed
                    continue
                ims, ext = res
                kargs['ext']=ext
                draw(ims, **kargs)
            return ims, ext
    ct=0
    stop=0
    for arg0 in args:
        if arg0 is None :
            continue
        if stop:
            break
        if isinstance(arg0, list):
            if len(arg0)==1:
                arg0=arg0[0]
        arg0=unwrap(arg0)
        if scipy.sparse.issparse(arg0):
            arg0=arg0.toarray()
        ct+=1
        if ct>1:
            plt.figure(figsize=kargs.get('figsize'))
        if isimg(arg0):
            #2-d numeric array
            img=np.squeeze(arg0)
            if img.size==0: #empty
                return
            elif len(img.shape) == 1 and img.shape[0]>0: #1-D array
                ncol = int(np.sqrt(img.shape[0]))
                nrow = int(img.shape[0] / ncol)
                if ncol*nrow==img.shape[0]:
                    img=img.reshape(ncol, nrow)
                else: #not a square image, plot as a curve
                    plt.plot(np.real(img))
                    return img
            try:
                plt.gca().images[-1].colorbar.remove()
            except (IndexError, AttributeError, KeyError, ValueError): #no previous image or colorbar
                pass
            #if img.dtype.name.find('complex')>-1:
            if np.iscomplexobj(img):
                img=np.real(img)
            if scipy.sparse.issparse(img):
                img=img.toarray()
            if 'center' in kargs:
                img=center(img, kargs['center'])
            if 'ext' not in kargs and 'dx' in kargs: #extent is [left, right, bottom, top]: x is along columns
                kargs['ext']=np.array([-img.shape[1]/2, img.shape[1]/2, -img.shape[0]/2, img.shape[0]/2])*kargs['dx']
            im=plt.imshow(img, extent=kargs.get('ext'), origin='lower', cmap='jet')
            #for im in gca().get_images():
            #im.set_clim(0, 0.5)
            if 'clim' in kargs and kargs['clim'] is not None:
                plt.clim(kargs['clim'])
            im_ratio = img.shape[0]/img.shape[1]
            if not 'colorbar' in kargs or kargs['colorbar']!=0:
                plt.colorbar(im, fraction=0.046*im_ratio, pad=0.04, extendfrac='auto')
            #cax = plt.gcf().add_axes([plt.gca().get_position().x1+0.01, plt.gca().get_position().y0, 0.02, plt.gca().get_position().height])
            #plt.colorbar(im, cax=cax)
            plt.grid(False)
            return img
            #if type(arg0) == list or arg0.dtype == object or arg0.ndim==3 or (isloc(arg0)==False and arg0.ndim==2 and arg0.shape[0]*100<arg0.shape[1]):  # list, array of array or 3d array
        else: #many images
            if isinstance(arg0, list) or isinstance(arg0, tuple):
                nframe = len(arg0)
            elif isinstance(arg0, np.ndarray):
                if iscell(arg0):
                    nframe = arg0.size
                    arg0=arg0.reshape(nframe)
                else:
                    if arg0.ndim==2 and arg0.shape[0]>arg0.shape[1]*100:
                        arg0=arg0.T #transpose
                    nframe = arg0.shape[0]
            else:
                nframe=0
                print('unknown type, skip', type(arg0))
            if nframe==0:
                continue

            if nframe>60:
                nframe=60
                print('Limit to first {} frames'.format(nframe))
            elif nframe==0:
                continue
            if 'nrow' in kargs:
                nrow=kargs['nrow']
                ncol=int(np.ceil(nframe/nrow))
            else:
                if 'nx' in kargs:
                    ncol=kargs['nx']
                elif 'ncol' in kargs:
                    ncol=kargs['ncol']
                elif nframe > 3:
                    ncol = int(np.ceil(np.sqrt(nframe)))
                else:
                    ncol = nframe
                if ncol>nframe:
                    ncol=nframe
                nrow = int(np.ceil(nframe/ncol))

            for iframe in range(nframe):
                if isimg(arg0[iframe]) and (ncol>1 or nrow > 1):
                    if iframe==0:
                        plt.clf()
                    plt.subplot(nrow, ncol, iframe+1)
                else:
                    plt.figure(figsize=kargs.get('figsize'))
                if isloc(arg0[iframe]) and len(args)>1:
                    draw(arg0[iframe], args[1][iframe], **kargs)
                    stop=1
                else:
                    draw(arg0[iframe], **kargs)
            plt.gcf().tight_layout() #fix layout
            return arg0
# Use as standalone script
if __name__ == "__main__":
    from readbin import readbin
    if len(sys.argv)<2:
        raise SystemExit('Usage: draw.py file1 [file2 ...]')
    args = list()
    for fn in sys.argv[1:]:
        args.append(readbin(fn))

    plt.figure()
    plt.clf()
    draw(*args)
    print('Close window to continue')
    plt.draw()
    plt.show()
