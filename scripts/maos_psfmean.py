#!/usr/bin/env python3
import os
import sys
import glob
import re
#import readbin
from astropy.io import fits
def maos_psfmean(fn1, fn2):
    fns=glob.glob(fn1)
    r=re.compile(fn2) 
    counts=dict();
    psfs=dict();
    header=dict();

    for fn in fns:
        tmp=r.split(fn)
        seed=tmp[1]
        if len(tmp)>3:
            thetax=tmp[2]
            thetay=tmp[3]
            isol=0
        else:
            thetax=0
            thetay=0
            isol=1
        try:
            hds=fits.open(fn)
        except OSError:
            print('unable to open ' + fn)
            continue
        with hds:
            nwvl=len(hds)
            key=(thetax, thetay)
            if key in counts:
                counts[key]+=1
                for iwvl in range(0, nwvl):
                    psfs[key][iwvl]+=hds[iwvl].data
            else:
                counts[key]=1
                psfs[key]=dict()
                header[key]=dict()
                for iwvl in range(0, nwvl):
                    psfs[key][iwvl]=hds[iwvl].data.astype(float) #copy, as the file is closed
                    header[key][iwvl]=hds[iwvl].header.copy()
    for key in psfs.keys():
        print(key, counts[key], ' seeds')
        hds=fits.HDUList()
        for iwvl in range(0, len(psfs[key])):
            psfs[key][iwvl]*=1./counts[key]
            #print("psfs sum to ", psfs[key][iwvl].sum())
            hdr=header[key][iwvl]
            if iwvl==0:
                hdr['EXTEND']=True
                hds.append(fits.PrimaryHDU(psfs[key][iwvl], hdr))
            else:
                hds.append(fits.ImageHDU(psfs[key][iwvl], hdr))
        if isol:
            fnout='evlpsfol_mean.fits'
        else:
            fnout='evlpsfcl_mean_x'+key[0]+'_y'+key[1]+'.fits'
        if os.path.isfile(fnout):
            os.remove(fnout)
        hds.writeto(fnout)
        
                
if __name__ == "__main__":
    pwd=os.getcwd()
    for arg in sys.argv[1:]:
        os.chdir(pwd)
        os.chdir(arg)
        print(arg)
        maos_psfmean("evlpsfcl_[0-9]*x*y*.fits", r'evlpsfcl_([0-9]+)_x([0-9.+-e]+)_y([0-9.+-e]+)\.fits')
        maos_psfmean("evlpsfol_[0-9]*.fits", r'evlpsfol_([0-9]+)\.fits')
        
