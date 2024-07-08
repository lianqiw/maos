#!/usr/bin/env python3
# Helper routines to parse MAOS results

import numpy as np
from numpy import nan
try:
    from natsort import natsorted
except:
    natsorted=sorted
import glob, os

from libaos import read
from readbin import readbin

def maos_cumu(files, seeds=None, nsim0=0): ##return cumulative average
    res,fds=maos_res(files,seeds,0,0)
    print(fds)
    nsim=res.shape[-1]
    if nsim0<=0:
        nsim0=min(5000, np.int(nsim*0.1))
    yy=np.arange(1, nsim+1-nsim0)
    xx=nsim0+yy
    yy.shape=(1,1,nsim-nsim0)
    resc=np.cumsum(res[:,:,nsim0:], axis=2)/yy
    resc[res[:,:,nsim0:]==0]=nan
    return resc,xx
def maos_res(fds, seeds=None, iframe1=0.2, iframe2=1, quiet=0):
    '''Results are in order of High, T/T, NGS total[, Focus]'''
    return maos_res_do(fds, "Res", seeds, iframe1, iframe2, quiet)
def maos_res_tot(fds, seeds=None, iframe1=0.2, iframe2=1, quiet=0):
    '''Results are High + Low tot'''
    res,fds=maos_res_do(fds, "Res", seeds, iframe1, iframe2, quiet)
    res2=res[:,0]+res[:,2]
    return res2,fds
def maos_res_hi(fds, seeds=None, iframe1=0.2, iframe2=1, quiet=0):
    '''Results are High order only'''
    res,fds=maos_res_do(fds, "Res", seeds, iframe1, iframe2, quiet)
    res2=res[:,0]
    return res2,fds
def maos_res_each(fds, seeds=None, iframe1=0.2, iframe2=1, quiet=0):
    '''Results fore all directions with dimension nfolder*ndirection*nmod. The modes are in PR, TT, PTTR.'''
    return maos_res_do(fds, "Resp", seeds, iframe1, iframe2, quiet)
def maos_res_do(fdin, name, seeds=None, iframe1=0.2, iframe2=1, quiet=0):
    '''Process maos results and average between firame1 to iframe2'''
    if type(fdin) is not list:
        fdin=[fdin]
    fds2=[]
    for fdini in fdin:
        fds2+=glob.glob(fdini+'/',recursive=1) #find only directory

    #natsorted not work well if there is trailing /
    fds2=natsorted([fd[0:-1] for fd in fds2])
    fds=[]
    resall=None
    nseed=0
    for fd in fds2: #loop over directory
        if seeds is None:
            fns=natsorted(glob.glob(fd+"/"+name+"_*.bin"))
        else:
            fns=list()
            for seed in seeds:
                fns.append(fd+'/{}_{}.bin'.format(name,seed))

        nseed=0
        mres=0
        split=-1
        for fn in fns: #loop over seed
            if not os.path.exists(fn):
                print(fn, 'does not exist')
                continue
            res=readbin(fn)
            if type(res) is tuple:
                res=res[0]
            if res is None or res.shape[0]==0:
                continue
            if name=="Res":
                if res.shape[0]>3 and res[3].size>0: #split tomography
                    res=res[3]
                    split=1
                else: #integrated
                    res=res[2]
                    res[:,0]=res[:,2] #convert from pr, tt, pttr to pttr, tt, tt to have the same as the split mode.
                    res[:,2]=res[:,1]
                    split=0
            elif name[-7:]=="Resclep": #per direction (old format)
                res=np.stack(res, axis=1)
                split=0
            elif name[-4:]=="Resp": #per direction (new format)
                res=np.stack(res[-1],axis=1) #last cell is clep
                split=0            
            else:
                print('Invalid result')
                continue
            n_valid=np.nonzero(res[:,0]>0)[0]
            if n_valid.size>0:
                n_valid=n_valid[-1]             
            else:
                continue
            if iframe1<1:
                n1=round(iframe1*n_valid)
            else:
                n1=iframe1

            if iframe2<=1:
                n2=round(iframe2*n_valid)
            else:
                n2=iframe2
            if n1 < n2 and n2 <= res.shape[0]:
                res=res[n1:n2]
                res=res.mean(0,keepdims=1)
            else:
                res=res.T
                res.shape=(1,res.shape[0],res.shape[1])
            #if max(res[-1]).all()==0:
            #    print(fn, 'Incomplete result')
            #    continue #unfinished result
    
            mres+=res*1e18
            nseed+=1
        if nseed>0:
            fds.append(fd)
            mres=mres*(1/nseed)
            if resall is None:
                resall=mres
            else:
                resall=np.vstack((resall, mres))
        elif quiet==0:
            print(fd, fns, nseed, ' has no valid results')
    if resall is None:
        resall=np.array([[nan,nan,nan]])
        #if not os.path.exists(fdin):
        #    print(fdin, 'does not exist')
        #else:
        if quiet==0:
            print(fdin, ' has no valid results')
    #if len(fds)>1:
    #    print(*fds, sep="\n")
    if quiet==0:
        print('{} results are read for {} seeds.'.format(len(fds), nseed))
    return (resall,np.array(fds))
    #else:
    #return resall
