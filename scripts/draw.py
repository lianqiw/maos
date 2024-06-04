#!/usr/bin/env python3
"""
This module contains routines to embed opds defined on coordinates to 2-d array and plot.
"""
#from aolib import *
import sys
import matplotlib.pyplot as plt
import numpy as np
import scipy
#from aotools import center #cyclic import
#import struct


def coord2grid(x, **kargs):
    xmin = x.min()
    xmax = x.max()
    x2 = x-xmin
    if x2.max() > 0:
        if 'dx' in kargs and kargs['dx']!=0:
            dx=kargs['dx']
        else:
            dx2=x2[x2>0].min() #assume regular
            dx=(xmax-xmin)*np.sqrt(np.pi*0.5/x.size) #assume evenly distributed
            if(abs(dx-dx2)<(dx+dx2)*0.1):
                dx=dx2
        dx2 = 1./dx
        nx = np.round((xmax-xmin)*dx2+1.).astype(int)
        ix = np.floor(x2*dx2+1e-5).astype(int)
    else:
        ix = np.zeros(x.shape)
        nx = 1
        dx = 1
    xi = [round(xmin/dx)*dx, round(xmax/dx)*dx]
    return (ix, nx, dx, xi)


def isvec(arg):
    return arg.ndim == 1 or (arg.ndim == 2 and (arg.shape(0) == 1 or arg.shape(1) == 1))


def isloc(arg):
    if type(arg) is not np.ndarray:
        return False
    if arg.dtype is object:
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

def locembed(loc, opd0, return_ext=0, **kargs):
    # draw(loc, opd): embed opd onto loc
    if loc.shape[0]>2 and loc.shape[1]==2: #need a transpose
        loc=loc.T
    opd=opd0.view()
    nloc = loc.shape[1]
    if opd.dtype==object: #cell
        ims = []
        ext = None
        for opdi in opd:
            if opdi.size>0:
                ims_i,ext=locembed(loc, opdi, return_ext, **kargs)
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
        im = np.full((nx*ny),np.NaN)
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
    if type(arg0) is scipy.sparse._csr.csr_matrix:
        return True
    elif type(arg0) is np.ndarray:
        while (arg0.dtype==object or arg0.ndim>2) and arg0.shape[0]==1:
            arg0=arg0[0]
        return arg0.dtype != object and (arg0.ndim==1 or (arg0.ndim==2 and arg0.shape[0]*100>arg0.shape[1] and arg0.shape[1]*100>arg0.shape[0]))
def draw(*args, **kargs):
    if(len(args)==0):
        return
    #if not 'keep' in kargs or kargs['keep'] == 0:
    #    plt.clf()
    arg0=args[0]
    if type(arg0) is np.ndarray:
        arg0=np.squeeze(arg0)
    if isloc(arg0):  # first argument is loc
        loc=arg0
        if loc.shape[0]>2 and loc.shape[1]==2: #need a transpose
            loc=loc.T
        if len(args) == 1:
            # plot grid
            plt.plot(loc[0], loc[1], '+')
            plt.axis('scaled')  # better than equal
            plt.xlabel('x (m)')
            plt.ylabel('y (m)')
        elif len(args) >= 2:
            ct=0
            for arg1 in args[1:]:
                ct+=1
                if ct>1:
                    plt.figure(figsize=kargs.get('figsize'))
                ims, ext2 = locembed(loc, arg1, return_ext=1, **kargs)
                kargs['ext']=ext2
                draw(ims, **kargs)
        return
    ct=0
    stop=0
    for arg0 in args:
        if arg0 is None :
            continue
        if stop:
            break
        if type(arg0) == list:
            if len(arg0)==1:
                arg0=arg0[0]
        if type(arg0) is np.ndarray:
            arg0=np.squeeze(arg0)
            while arg0.shape[0]==1 and (arg0.dtype == object or arg0.ndim>2):
                arg0=np.squeeze(arg0[0])
        if type(arg0) is scipy.sparse._csr.csr_matrix:
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
                else:
                    raise(Exception('Unable to reshape 1d array'))
            try:
                plt.gca().images[-1].colorbar.remove()
            except:
                pass
            #if img.dtype.name.find('complex')>-1:
            if np.iscomplexobj(img):
                img=np.real(img)
            if type(img) is scipy.sparse._csr.csr_matrix:
                img=img.toarray()
            if 'center' in kargs:
                img=center(img, kargs['center'])
            if 'ext' not in kargs and 'dx' in kargs:
                kargs['ext']=np.array([-img.shape[0]/2, img.shape[0]/2, -img.shape[1]/2, img.shape[1]/2])*kargs['dx']
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
            #if type(arg0) == list or arg0.dtype == object or arg0.ndim==3 or (isloc(arg0)==False and arg0.ndim==2 and arg0.shape[0]*100<arg0.shape[1]):  # list, array of array or 3d array
        else: #many images
            if type(arg0) == list or type(arg0)==tuple:
                nframe = len(arg0)
            elif type(arg0) == np.ndarray:
                if arg0.dtype==object:
                    nframe = arg0.size
                    arg0=arg0.reshape(nframe)
                else:
                    if arg0.ndim==2 and arg0.shape[0]>arg0.shape[1]*100:
                        arg0=arg0.T #transpose
                    nframe = arg0.shape[0]
            else:
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

# Use as standalone script
if __name__ == "__main__":
    from readbin import read
    args = list()
    for fn in sys.argv[1:]:
        args.append(read(fn))

    plt.figure()
    plt.clf()
    draw(*args)
    print('Close window to continue')
    plt.draw()
    plt.show()
