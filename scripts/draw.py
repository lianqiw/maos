#!/usr/bin/env python3
"""
This module contains routines to embed opds defined on coordinates to 2-d array and plot.
"""
#from aolib import *
import matplotlib.pyplot as plt
import numpy as np
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
    return arg.ndim == 2 and ((arg.shape[0]>2 and arg.shape[1]==2) or (arg.shape[0] == 2 and arg.shape[1]>2))


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
    opd.shape=(nframe, nloc) #reshape() may do copy. assign .shape prevents copy

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

def draw(*args, **kargs):
    #if not 'keep' in kargs or kargs['keep'] == 0:
    #    plt.clf()
    arg0=args[0]
    if arg0 is None:
        return
    if type(arg0) == list:
        if len(arg0)==1:
            arg0=arg0[0]
    elif arg0.dtype == object or arg0.ndim==3:
        if arg0.shape[0]==1:
            arg0=arg0[0]
    if type(arg0) == list or arg0.dtype == object or arg0.ndim==3:  # list, array of array or 3d array
        kargs['keep'] = 1  # do not clear
        if type(arg0) == list:
            nframe = len(arg0)
        elif type(arg0) == np.ndarray and arg0.dtype==object:
            old_shape=arg0.shape
            nframe = arg0.size
            arg0.shape=(arg0.size,)
        elif arg0.ndim==3:
            nframe = arg0.shape[0]
        else:
            raise(Exception('Unknow data type'))
        if nframe>60:
            nframe=60
            print('Limit to first {} frames'.format(nframe))
        elif nframe==0:
            return
        if 'nx' in kargs:
            nx=kargs['nx']
        elif nframe > 3:
            nx = int(np.ceil(np.sqrt(nframe)))
        else:
            nx = nframe
        if nx>nframe:
            nx=nframe
        ny = int(np.ceil(nframe/nx))
        # print(nx,ny)
        if nx>1 or ny > 1:
            plt.clf()
        for iframe in range(nframe):
            if nx>1 or ny > 1:
                plt.subplot(ny, nx, iframe+1)
            if len(args) == 1:
                draw(arg0[iframe], **kargs)
            elif len(args) == 2:
                #if args[1].dtype==object or args[1].ndim>1:
                draw(arg0[iframe], args[1][iframe],  **kargs)
                #else:
                #    draw(arg0[iframe], args[1],  **kargs)
            else:
                print('Too many arguments')
        if type(arg0)==object:
            arg0.shape=old_shape
    elif isloc(arg0):  # first argument is loc
        if arg0.shape[0]>2 and arg0.shape[1]==2: #need a transpose
            arg0=arg0.T
        if len(args) == 1:
            # plot grid
            plt.plot(arg0[0], arg0[1], '+')
            plt.axis('scaled')  # better than equal
            plt.xlabel('x (m)')
            plt.ylabel('y (m)')
        elif len(args) == 2:
            # draw(loc, opd): embed opd onto loc
            ims, ext2 = locembed(arg0, args[1], 1, **kargs)
            # print('ext2=',ext2)
            kargs['ext']=ext2
            draw(ims, **kargs)
        else:
            print('Too many arguments')
    
    else: #2-d numeric array
        img = np.squeeze(arg0)
        if img.size==0: #empty
            return
        elif len(img.shape) == 1 and img.shape[0]>0: #1-D array
            nx = int(np.sqrt(img.shape[0]))
            ny = int(img.shape[0] / nx)
            if nx*ny==img.shape[0]:
                img.shape=(nx, ny)
            else:
                raise(Exception('Unable to reshape 1d array'))
        try:
            plt.gca().images[-1].colorbar.remove()
        except:
            pass
        if img.dtype.name.find('complex')>-1:
            img=np.real(img)
        im=plt.imshow(img, extent=kargs.get('ext'), origin='lower', cmap='jet')
        #for im in gca().get_images():
        #im.set_clim(0, 0.5)
        if 'clim' in kargs and kargs['clim'] is not None:
            plt.clim(kargs['clim'])
        im_ratio = img.shape[0]/img.shape[1]
        plt.colorbar(im, fraction=0.046*im_ratio, pad=0.04, extendfrac='auto')
        #cax = plt.gcf().add_axes([plt.gca().get_position().x1+0.01, plt.gca().get_position().y0, 0.02, plt.gca().get_position().height])
        #plt.colorbar(im, cax=cax)
        plt.grid(False)

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
