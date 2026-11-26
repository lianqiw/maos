#!/usr/bin/env python3
"""Average the PSFs saved by maos over seeds.

maos saves the time averaged PSF of each seed in evlpsfcl_<seed>_x<thetax>_y<thetay>.fits
(closed loop, per evaluation direction) and evlpsfol_<seed>.fits (open loop),
with one image (HDU) per wavelength. This script averages them over all the
seeds found and writes evlpsfcl_mean_x<thetax>_y<thetay>.fits and
evlpsfol_mean.fits in the same folder, keeping the header of the first seed.
Existing output files are overwritten.

Usage:
    maos_psfmean.py folder1 [folder2 ...]
    or in python: maos_psfmean(glob_pattern, regex) in the current directory.

Requirements: astropy.
"""
import os
import sys
import glob
import re
#import readbin
from astropy.io import fits
def maos_psfmean(fn1, fn2):
    """Average the PSF fits files matching fn1 over seeds.

    Args:
        fn1: glob pattern of the input files.
        fn2: regular expression that matches the whole file name with groups
            (seed) for open loop or (seed, thetax, thetay) for closed loop.
            Files that do not match are skipped.

    Writes evlpsfol_mean.fits or evlpsfcl_mean_x<thetax>_y<thetay>.fits in
    the current directory. Files that cannot be read, or whose number of
    images differs from the first file of the same direction, are skipped.
    """
    fns=glob.glob(fn1)
    r=re.compile(fn2) 
    counts=dict();
    psfs=dict();
    header=dict();

    isol=1
    for fn in fns:
        tmp=r.split(fn)
        if len(tmp)<3 or r.fullmatch(fn) is None: #not matching the expected file name
            print('skip ' + fn)
            continue
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
            if any(hd.data is None for hd in hds): #empty image
                print('skip {}: has empty image'.format(fn))
                continue
            if key in counts and nwvl!=len(psfs[key]):
                print('skip {}: has {} images, expect {}'.format(fn, nwvl, len(psfs[key])))
                continue
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
    if len(sys.argv)<2:
        raise SystemExit('Usage: maos_psfmean.py folder1 [folder2 ...]')
    pwd=os.getcwd()
    for arg in sys.argv[1:]:
        os.chdir(pwd)
        os.chdir(arg)
        print(arg)
        #thetax and thetay are numbers such as 0, -17, 1.5, or 1e-05
        maos_psfmean("evlpsfcl_[0-9]*x*y*.fits", r'evlpsfcl_([0-9]+)_x([-0-9.+e]+)_y([-0-9.+e]+)\.fits')
        maos_psfmean("evlpsfol_[0-9]*.fits", r'evlpsfol_([0-9]+)\.fits')
        
