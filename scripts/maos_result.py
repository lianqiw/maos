#!/usr/bin/env python3
"""Helper routines to parse MAOS results.

Read the wavefront error results (Res_<seed>.bin, Resp_<seed>.bin) saved by
maos in one or more result folders, average them over a range of time steps
and over the seeds, and return the mean wavefront variance in nm^2.

Usage:
    from maos_result import *
    res,fds=maos_res(fds, seeds=None, iframe1=0.2, iframe2=1, quiet=1)      #field averaged WFE
    res,fds=maos_res_tot(fds, seeds=None, iframe1=0.2, iframe2=1, quiet=1)  #high order + low order total
    res,fds=maos_res_hi(fds, seeds=None, iframe1=0.2, iframe2=1, quiet=1)   #high order only
    res,fds=maos_res_each(fds, seeds=None, iframe1=0.2, iframe2=1, quiet=1) #each evaluation direction
    resc,steps=maos_cumu(fds, seeds=None, nsim0=0)                          #cumulative average

Common arguments:
    fds: a folder name or a list of folder names. Wildcards are expanded with glob.
    seeds: list of seeds. Default: all the result files in each folder.
    iframe1, iframe2: time steps to average. A value <1 (<=1 for iframe2) is a
        fraction of the number of completed steps; otherwise a step index.
        With iframe1=0 and iframe2=0 the time history is returned instead.
    quiet: 0 to print diagnostic messages.

Returns of maos_res* functions:
    (res, fds): res has one row per folder with valid results (a single row of
    NaN if there is none); fds is a numpy array of the folders used.

Requirements: numpy, readbin.py in the module search path; natsort is optional.
See also docs/40_results.md.
"""

import numpy as np
from numpy import nan
try:
    from natsort import natsorted
except ImportError:
    natsorted=sorted
import glob, os

from readbin import readbin

def maos_cumu(files, seeds=None, nsim0=0): ##return cumulative average
    """Cumulative average of the field averaged WFE (see maos_res) versus time step.

    Args:
        files: folder name(s), see maos_res.
        seeds: list of seeds (default: all).
        nsim0: number of initial steps to skip. Default (<=0): 10% of the
            steps, up to 5000.

    Returns:
        (resc, xx): resc of shape (nfolder, ncol, nsim-nsim0) (ncol as in maos_res) is the cumulative
        average in nm^2 from step nsim0 (NaN where the result is 0, i.e., not
        yet computed); xx is the corresponding number of steps (end step index).
    """
    res,fds=maos_res(files,seeds,0,0)
    print(fds)
    nsim=res.shape[-1]
    if nsim0<=0:
        nsim0=min(5000, int(nsim*0.1))
    yy=np.arange(1, nsim+1-nsim0)
    xx=nsim0+yy
    yy.shape=(1,1,nsim-nsim0)
    resc=np.cumsum(res[:,:,nsim0:], axis=2)/yy
    resc[res[:,:,nsim0:]==0]=nan
    return resc,xx
def maos_res(fds, seeds=None, iframe1=0.2, iframe2=1, quiet=1):
    '''Field averaged closed loop WFE (nm^2) of shape (nfolder, 3[+1]).

    Results are in order of High, T/T, NGS total[, Focus]. For integrated
    tomography, High is the piston/tip/tilt removed WFE and both T/T and NGS
    total are the tip/tilt WFE. Focus is included when the simulation has a
    separate focus mode. With iframe1=iframe2=0, the time history of
    shape (nfolder, ncol, nstep) is returned.
    '''
    return maos_res_do(fds, "Res", seeds, iframe1, iframe2, quiet)
def maos_res_tot(fds, seeds=None, iframe1=0.2, iframe2=1, quiet=1):
    '''Results are High + Low tot: total WFE (nm^2) of shape (nfolder,). See maos_res.'''
    res,fds=maos_res_do(fds, "Res", seeds, iframe1, iframe2, quiet)
    res2=res[:,0]+res[:,2]
    return res2,fds
def maos_res_hi(fds, seeds=None, iframe1=0.2, iframe2=1, quiet=1):
    '''Results are High order only: WFE (nm^2) of shape (nfolder,). See maos_res.'''
    res,fds=maos_res_do(fds, "Res", seeds, iframe1, iframe2, quiet)
    res2=res[:,0]
    return res2,fds
def maos_res_each(fds, seeds=None, iframe1=0.2, iframe2=1, quiet=1):
    '''Results for all directions with dimension nfolder*ndirection*nmod (nm^2).

    The modes are in PR, TT, PTTR (closed loop piston removed, tip/tilt, and
    piston/tip/tilt removed WFE, from the last part of Resp_<seed>.bin).
    '''
    return maos_res_do(fds, "Resp", seeds, iframe1, iframe2, quiet)
def maos_res_do(fdin, name, seeds=None, iframe1=0.2, iframe2=1, quiet=1):
    '''Process maos results and average between iframe1 to iframe2.

    Args:
        fdin: folder name(s) with glob wildcards.
        name: "Res" (field average) or "Resp" (each direction).
        seeds, iframe1, iframe2, quiet: see the module docstring.

    Returns:
        (res, fds): averaged results in nm^2 and the folders used. See
        maos_res and maos_res_each for the shape of res.
    '''
    if not isinstance(fdin, (list, tuple)):
        fdin=[fdin]
    fds2=[]
    for fdini in fdin:
        fds2+=glob.glob(os.fspath(fdini)+'/',recursive=1) #find only directory

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
            if isinstance(res, tuple):
                res=res[0]
            if res is None or res.shape[0]==0:
                continue
            #determine result cell
            #Res has 4 cells: OL, None, CL, AHST
            if name=="Res":
                if res.shape[0]>3 and res[3].size>0: #split tomography
                    res=res[3] #AHST: LGS, TT, NGS(TT+PS)
                    split=1
                else: #integrated
                    res=res[2] #PR, TT, PTTR.
                    res[:,0]=res[:,2] #Conver to AHST format
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
            #determine valid simulation duration in steps: up to the last step with a result
            n_valid=np.nonzero(res[:,0]>0)[0]
            if n_valid.size>0:
                n_valid=n_valid[-1]+1
            else:
                continue
            if iframe1<1: #fraction of valid steps
                n1=round(iframe1*n_valid)
            else: #step index
                n1=int(iframe1)

            if iframe2<=1:
                n2=round(iframe2*n_valid)
            else:
                n2=int(iframe2)
            if n1==0 and n2==0: #want time history
                res=res.T
                res.shape=(1,res.shape[0],res.shape[1])
            elif n1 < n2 and n2 <= res.shape[0]:
                res=res[n1:n2]
                res=res.mean(0,keepdims=1)
            else:
                print('Invalid configuration for {}: n1={},n2={},n_valid={}'.format(fn, n1, n2, n_valid))
                res=np.full((1,)+res.shape[1:], np.nan) #same shape as the averaged result
            #if max(res[-1]).all()==0:
            #    print(fn, 'Incomplete result')
            #    continue #unfinished result
            mres+=res*1e18 #m^2 to nm^2
            nseed+=1
        if nseed>0:
            fds.append(fd)
            mres=mres*(1/nseed)
            if resall is None:
                resall=mres
            elif resall.shape[1:]!=mres.shape[1:]:
                raise ValueError(f'{fd}: result of shape {mres.shape[1:]} does not match that of the previous folders {resall.shape[1:]}. Process them separately.')
            else:
                resall=np.vstack((resall, mres))
        elif quiet==0:
            print(fd, fns, nseed, ' has no valid results')
    if resall is None:
        resall=np.array([[nan,nan,nan]])
        if name!="Res":
            resall.shape=(1,1,3)
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
