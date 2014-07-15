import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as pyplot
import math
def draw_each(image, iscell):
    shape=image.shape
    if len(shape)==2 and min(shape)==1:
        shape=(max(shape)) #vector
    if len(shape)==2 and min(shape)>2:
        pyplot.imshow(imi.astype(np.double), origin='lower', interpolation='nearest')
        if iscell:
            pyplot.axis('off')
        else:
            pyplot.colorbar()
    elif len(shape)==2 and shape[0]=='2':
        pyplot.plot(image[0], image[1], '+')
        
def draw(image, ifig=0):
    fig=pyplot.figure(ifig, figsize=(8,6), dpi=120)
    iscell=0
    nrow=1
    ncol=1
    while image.dtype==np.dtype(object) and len(image.flat)==1:
        image=image.flat[0]
    if image.dtype==np.dtype(object): #cell array
        iscell=1
        shape=image.shape
        if len(shape)==2:
            nrow=shape(1)
            ncol=shape(2)
        else:
            nrow=1
            ncol=1
            for i in shape:
                ncol=ncol*i
        if nrow == 1 or ncol == 1:
            tot=nrow*ncol
            if tot>10:
                ncol=math.sqrt(tot)
                if ncol>6:
                    ncol=6
                nrow=tot/ncol
        if ncol>10:
            ncol=10
        if nrow>4:
            nrow=4
    for ip in xrange(0, nrow*ncol):
        spl=fig.add_subplot(nrow, ncol, ip)
        if iscell:
            draw_each(image.flat[ip], iscell)
        else:
            draw_each(image, iscell)
    
    return fig
