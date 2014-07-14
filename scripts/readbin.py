#/usr/bin/env python3
from numpy import *
import gzip
import struct
isfits=0


def readuint16(fp):
    return struct.unpack("<H", fp.read(2))[0]
def readbin(file):
    if file[-5:]=='.fits' or file[-8:] == '.fits.gz':
        isfits=1
    elif file[-4:]=='.bin' or file[-7:] == '.bin.gz':
        isfits=0
    else:
        print 'Unknown file name, assume .bin'
        isfits=0
    try:
        fp=open(file, 'rb')
        magic=fp.read(2);
        if magic==0x8b1f:
            fp.close()
            fp=gzip.open(file,'rb')
        else:
            fp.seek(0, 0)
    
        err=0
        count=0
        out=list()
        header=list()
        if isfits:
            while err==0:
                [out_i, header_i, err]=readbin_do(fp, isfits)
                if err==0:
                    out.append(out_i)
                    header.append(header_i)
        else:
            [out, header, err]=readbin_do(fp, isfits)
    finally:
        fp.close()

    return (out, header)

def readbin_do(fp, isfits):
    if isfits:
        [magic, nx, ny, header]=readfits_header(fp)
    else:
        [magic, nx, ny, header]=readbin_header(fp)

bitpix2magic={-32:0x6408, -64:0x6402, 32:0x6405, 16:0x640B, 8:0x640, 0:0}
def readfits_header(fp):
    END=0
    page=0
    bitpix=0
    nx=0
    ny=0
    header=''
    while !END:
        if page==0:
            try:
                res=fp.read(80)
            except:
                break
            if res[0:6] !='SIMPLE' and res[0:16] !='XTENSION= ''IMAGE':
                print('Unknown data format in filts file');
                break
            res=fp.read(80)
            bitpix=int(res[10:30])
            res=fp.read(80)
            naxis=int(res[10:30])
            if naxis>2:
                print('Data type not supported')
                break
            if naxis>0:
                res=fp.read(80)
                nx=int(res[10:30])
            if naxis>1:
                res=fp.read(80)
                ny=int(res[10:30])
            for i in range(3+naxis, 36):
                res=fp.read(80)
                if res[0:7]=='COMMENT':
                    header+=res[10:].rstrip()
                elif res[0:3]=='END':
                    END=1
            page=page+1
    
