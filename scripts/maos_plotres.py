#!/bin/env ipython3
try:
    from startup import *
except:
    pass
import sys
import numpy as np
import matplotlib.pyplot as plt
from readbin import readbin
import glob

def maos_plotres(fns):
    legs=list()
    res=list()
    print("\n".join(fns))

    for fn in fns:
        fn2=glob.glob(fn+"/Res_*.bin")
        for fn3 in fn2:
            ires=readbin(fn3)[3]
            ires[:,2]=ires[:,2]-ires[:,1]-ires[:,3];
            ires=np.sqrt(ires)*1e9;
            res.append(ires)
            legs.append(fn)
    ylabel=['High','TT','PS','Focus']
    for ic in range(0,4):
        plt.figure(ic)
        plt.clf()
        for ires in res:
            plt.plot(ires[:,ic])
        plt.legend(legs)
        plt.xlabel('Time step')
        plt.ylabel(ylabel[ic]+ " WFE (nm)")
    plt.show(block=True)
    return (res,legs)
if __name__ == "__main__":
    fns=list()
    for dirs in sys.argv[1:]:
        fns+=glob.glob(dirs)
    (res,legs)=maos_plotres(fns)
