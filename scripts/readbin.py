#!/usr/bin/env python3
import sys
import numpy as np
import scipy.sparse as sparse
import os
import gzip
import struct
import socket
# magic number are used in bin files to indicate data type
# dname descrinebs data type. M_ types are fundamental types. 
# MC_ types has been superseded by MCC_ANY
magic2dname={
    25600: 'M_CSP64',
    25601: 'M_SP64',
    25602: 'M_DBL',
    25603: 'M_INT64',
    25604: 'M_CMP',
    25605: 'M_INT32',
    25606: 'M_CSP32',
    25607: 'M_SP32',
    25608: 'M_FLT',
    25609: 'M_ZMP',
    25610: 'M_INT8',
    25611: 'M_INT16',
    25616: 'MC_CSP',
    25617: 'MC_SP',
    25618: 'MC_DBL',
    25619: 'MC_INT64',
    25620: 'MC_CMP',
    25621: 'MC_INT32',
    25633: 'MCC_ANY',
    25634: 'MCC_DBL',
    25636: 'MCC_CMP',
    #25856: 'M_COMMENT',
    #26112: 'M_SKIP'
}

dname2type={
    'M_ZMP': np.complex64,
    'M_CMP': np.complex128,
    'M_CSP32': np.complex64,
    'M_CSP64': np.complex128,
    'M_DBL': np.double,
    'M_FLT': np.float32,
    #'M_COMMENT': object,
    'M_INT16': np.int16,
    'M_INT32': np.int32,
    'M_INT64': np.int64,
    'M_INT8': np.int8,
    'M_SKIP': object,
    'M_SP32': np.float32,
    'M_SP64': np.double
}

bitpix2magic={
    -32:0x6408,
    -64:0x6402,
    32:0x6405,
    16:0x640B,
    8:0x640,
    0:0
}

def readuint16(fp):
    return struct.unpack("<H", fp.read(2))[0]
def readuint32(fp):
    return struct.unpack("<I", fp.read(4))[0]
def readuint64(fp):
    return struct.unpack("<Q", fp.read(8))[0]
def readvec(fp, datatype, nxy):#enable read from socket and file
    if fp.seekable():
        return np.fromfile(fp, dtype=datatype, count=nxy, sep='')
    else:#for socket reading, fromfile fails.
        buf=fp.read(nxy*datatype.itemsize)
        return np.frombuffer(buf, dtype=datatype)

def readbin(file, want_header=0):
    isfits=False
    isfile=True
    out=np.array(())
    header=b''
    if isinstance(file, socket.socket):
        file=file.fileno()
        isfits=False
        isfile=False
    elif not os.path.isfile(file): #file not found
        if os.path.isfile(file+'.bin'):
            file=file+'.bin'
        elif os.path.isfile(file+'.bin.gz'):
            file=file+'.bin.gz'
        elif os.path.isfile(file+'.fits'):
            file=file+'.fits'
            isfits=True
        elif os.path.isfile(file+'.fits.gz'):
            file=file+'.fits.gz'
            isfits=True
        else:
            print('File does not exist:', file)
    if os.path.isfile(file) or not isfile:
        if (file[-5:]=='.fits' or file[-8:] == '.fits.gz'):
            isfits=True
        try:
            fp=open(file, 'rb', closefd=isfile)
            if isfile:
                magic=readuint16(fp)
                if magic==0x8b1f:
                    fp.close()
                    fp=gzip.open(file,'rb')
                else:
                    fp.seek(0, 0)
            (out, header, err)=readbin_auto(fp, isfits)
        except Exception as error:
            print("readbin failed:", file, error)
        finally:
            fp.close()
    if want_header: 
        return (out, header)
    else:
        return out
def readbin_auto(fp, isfits, scan=0):
    err=0
    out=list()
    header=list()
    if scan or isfits:
        while err==0:
            [out_i, header_i, err]=readbin_do(fp, isfits)
            if err==0:
                out.append(out_i)
                header.append(header_i)
        if len(out)==1:
            out=out[0]
            header=header[0]
        else:
            out=np.array(out)
    else:
        [out, header, err]=readbin_do(fp, isfits)
    if err<0:
        print('readbin: error happened', err)
    return (out, header, err)

def readbin_do(fp, isfits):
    err=0
    header=''
    if isfits:
        [magic, nx, ny, header]=readfits_header(fp)
    else:
        [magic, nx, ny, header]=readbin_header(fp)

    if magic<0:
        dname='Unknown'
        err=-1
    elif magic==0:
        dname='EOF'
        err=1
    else:
        try:
            dname=magic2dname[magic & 0xFFFF]
        except:
            raise ValueError("Invalid magic number {}".format(magic))
        
    if dname[0:2]=='MC': #cell array
        if nx>0 and ny>0:
            header_cell=header
            header=np.zeros((ny, nx), dtype=object)
        else:
            return readbin_auto(fp, isfits, 1)
        out=np.zeros((ny, nx), dtype=object)

        for iy in range(0, ny):
            for ix in range(0, nx):
                out[iy,ix], header[iy,ix], err=readbin_do(fp, isfits)
                if len(header[iy,ix])==0:
                    header[iy,ix]=header_cell
                if err:
                    break
            if err:
                break
    
        if ny==1:
            out=out[0,]
            header=header[0,]
    elif 'SP' in dname and nx>0 and ny>0: #sparse matrix
        datatype=np.dtype(dname2type[dname])
        nz=readvec(fp, np.dtype(np.int64), 1)[0]
        #nz=np.fromfile(fp, dtype=np.uint64, count=1, sep='')[0]
        if '64' in dname:
            ijtype=np.dtype(np.int64)
        else:
            ijtype=np.dtype(np.int32)
        Jc=readvec(fp, ijtype, ny+1).astype(int)
        Ir=readvec(fp, ijtype, nz).astype(int)
        P=readvec(fp, datatype, nz)
        out=sparse.csr_matrix((P, Ir, Jc), shape=(ny, nx))
    elif dname[0:2]=='M_' and nx>0 and ny>0: #dense matrix
        datatype=np.dtype(dname2type[dname])
        if isfits:
            datatype=datatype.newbyteorder('>')
        out=readvec(fp, datatype, nx*ny)
        if len(out) != nx*ny:
            print('Wrong length is read. Got', len(out), ' expect ', nx*ny, fp.name)
        elif ny>1:
            out.shape=(ny, nx)
        
        if isfits:
            byteread=datatype.itemsize*nx*ny
            byteleft=byteread%2880
            if byteleft:
                fp.seek(2880-byteleft,1) #skip
    else:
        out=np.array(())
    return (out, header, err)


def readfits_header(fp):
    END=0
    page=0
    bitpix=0
    nx=0
    ny=0
    header=b''
    while END==0:
        start=0
        if page==0: #first page
            res=fp.read(80)
            if len(res)!=80: #end of line
                END=2
                break
            elif res[0:6] !=b'SIMPLE' and res[0:16] != b"XTENSION= 'IMAGE":
                END=3
                print('Unknown data format in fits file', END, page, start, len(res), fp.tell(), res)
                break
            res=fp.read(80)
            bitpix=int(res[10:30])
            res=fp.read(80)
            naxis=int(res[10:30])
            if naxis>2:
                print('Data type not supported')
                END=4
                break
            if naxis>0:
                res=fp.read(80)
                nx=int(res[10:30])
            if naxis>1:
                res=fp.read(80)
                ny=int(res[10:30])
            start=3+naxis
        for i in range(start, 36):
            res=fp.read(80)
            if res[0:7]==b'COMMENT':
                header+=res[10:].rstrip()
            elif res[0:3]==b'END':
                END=1 #do not break. let all data be read
        page=page+1
    if END == 1:
        return (bitpix2magic[bitpix], nx, ny, header)
    elif END == 2: #end of file
        return (0, 0, 0, '')
    else:
        return (-1, 0, 0, '')
def readbin_magic(fp):
    magic=readuint32(fp)
    M_SKIP=26112
    if magic==M_SKIP: #padding
        magic=readuint32(fp)
    return magic
def readbin_header(fp):
    M_COMMENT=25856
    header=b''
  
    magic=readbin_magic(fp)
    while magic&0xFFFF==M_COMMENT:
        nlen=readuint64(fp)
        header+=fp.read(nlen).strip() #.decode('utf-8')
        nlen2=readuint64(fp)
        magic2=readbin_magic(fp)
        if nlen!=nlen2 or magic!=magic2:
            raise NameError('Header verification failed')
        magic=readbin_magic(fp)
    nx=readuint64(fp)
    ny=readuint64(fp)
   
    return (magic, nx, ny, header)

if __name__ == '__main__':
    if len(sys.argv)>1:
        for fn in sys.argv[1:]:
            res=readbin(fn)
            print(fn, ' is ', res.shape)
    else:
        raise ValueError("Usage: res=readbin.readbin(file name)")
