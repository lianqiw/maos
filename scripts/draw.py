#!/usr/bin/env python3
"""
This module contains routines to embed opds defined on coordinates to 2-d array and plot.
"""
#from aolib import *
import matplotlib.pyplot as plt
import numpy as np
#import struct


def coord2grid(x):
    xmin = x.min()
    xmax = x.max()
    x2 = x-xmin
    if x2.max() > 0:
        dx = x2[x2 > 0].min()
        dx2 = 1/dx
        nx = np.round((xmax-xmin)*dx2+1.).astype(int)
        ix = np.floor(x2*dx2+1e-5).astype(int)
    else:
        ix = np.zeros(x.shape)
        nx = 1
        dx = 1
    xi = [xmin, xmax]
    return (ix, nx, dx, xi)


def isvec(arg):
    return arg.ndim == 1 or (arg.ndim == 2 and (arg.shape(0) == 1 or arg.shape(1) == 1))


def isloc(arg):
    return arg.ndim == 2 and arg.shape[0] == 2


def locembed(loc, opd0):
    # draw(loc, opd): embed opd onto loc
    opd=opd0.view()
    (ix, nx, dx, xi) = coord2grid(loc[0])
    (iy, ny, dy, yi) = coord2grid(loc[1])

    nloc = loc.shape[1]
    if len(opd.shape)==1: #vector
        nframe = opd.size/nloc
        opd.shape=(int(nframe), nloc) #reshape() may do copy. assign .shape prevents copy
    elif opd.shape[0]==nloc and opd.shape[1]!=nloc:
        opd=opd.T
    elif opd.shape[0]!=nloc and opd.shape[1]==nloc:
        pass
    else:
        print("data has wrong shape", opd.shape)
        return None,None
    nframe=opd.shape[0]
    
    # print('opd=',opd.shape)
    ims = np.empty(nframe, dtype=object)
    for iframe in range(nframe):
        im = np.full((nx*ny),np.NaN)
        im[ix+iy*nx] = opd[iframe, :]
        im.shape = (ny, nx)  # fixed 2020-10-09
        ims[iframe] = im
    ext = np.r_[xi[0]-dx/2, xi[1]+dx/2, yi[0]-dy/2, yi[1]+dy/2]
    if nframe > 1:
        return ims, ext
    else:
        return im, ext


def draw(*args, **kargs):
    #if not 'keep' in kargs or kargs['keep'] == 0:
    #    plt.clf()
    if type(args[0]) == None:
        return
    if type(args[0]) == list or args[0].dtype == object:  # array of array
        kargs['keep'] = 1  # do not clear
        if type(args[0]) == list:
            nframe = len(args[0])
        else:
            nframe = args[0].size
        if nframe > 3:
            nx = np.ceil(np.sqrt(nframe))
        else:
            nx = nframe
        ny = int(np.ceil(nframe/nx))
        # print(nx,ny)
        for iframe in range(nframe):
            if nx*ny > 1:
                plt.subplot(ny, nx, iframe+1)
            if len(args) == 1:
                draw(args[0][iframe], **kargs)
            elif len(args) == 2:
                draw(args[0][iframe], args[1][iframe],  **kargs)
            else:
                print('Too many arguments')
    elif isloc(args[0]):  # first argument is loc
        if len(args) == 1:
            # plot grid
            plt.plot(args[0][0], args[0][1], '+')
            plt.axis('scaled')  # better than equal
            plt.xlabel('x (m)')
            plt.ylabel('y (m)')
        elif len(args) == 2:
            # draw(loc, opd): embed opd onto loc
            ims, ext2 = locembed(args[0], args[1])
            # print('ext2=',ext2)
            draw(ims, ext=ext2)
            return ims
        else:
            print('Too many arguments')
    elif args[0].ndim > 2:
        kargs['keep'] = 1  # do not clear
        nframe = args[0].shape[0]
        if nframe > 3:
            nx = np.ceil(np.sqrt(nframe))
        else:
            nx = nframe
        ny = int(np.ceil(nframe/nx))
        # print(nx,ny)
        for iframe in range(nframe):
            if nx*ny > 1:
                plt.subplot(ny, nx, iframe+1)
            draw(args[0][iframe, ], **kargs)
    else:
        img = np.squeeze(args[0])
        if len(img.shape) == 1:
            nx = int(np.sqrt(img.shape[0]))
            ny = int(img.shape[0] / nx)
            if nx*ny==img.shape[0]:
                img.shape=(nx, ny)
            else:
                raise(Exception('Unable to reshape 1d array'))
        if 'ext' in kargs:
            plt.imshow(img, extent=kargs['ext'], origin='lower', cmap='jet')
        else:
            plt.imshow(img, origin='lower', cmap='jet')
        plt.colorbar()
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
