#!/usr/bin/env python3
"""Pure Python reader for the MAOS .bin and .fits file formats.

MAOS saves numerical arrays, sparse matrices, and (nested) cell arrays in its
own .bin format (see sys/bin.c and sys/bin.h) or as FITS images. This module
reads them into numpy/scipy objects without needing the compiled aolib.so
(see interface.py/libaos.py for the ctypes based alternative). writebin.py
writes them back.

Usage:
    from readbin import readbin
    res=readbin('Res_1')  #tries Res_1, Res_1.bin, Res_1.bin.gz, Res_1.fits, Res_1.fits.gz
    python3 readbin.py file1.bin [file2.fits ...]  #print the shape of each file

Conventions:
    - MAOS arrays are column major. An array of nx*ny in MAOS is returned with
      shape (ny, nx); one with ny==1 is returned as a 1-d array of shape (nx,).
    - A cell array is returned as a numpy object array with the same convention.
      Empty elements are returned as empty arrays.
    - Sparse matrices (MAOS CSC) are returned as scipy.sparse.csr_matrix of
      shape (ny, nx), i.e., the transpose convention as dense arrays.
    - Each returned dense or sparse array carries the header in its ``header``
      attribute: a string (or numpy array of strings, one per line) for .bin
      files, a dict of keyword/value strings for .fits files.
    - A FITS file with several images (HDUs) is returned as an array of them.

Supported compression: gzip (detected from the content, regardless of the file
name). zstd compressed files (.zst), which MAOS can also write, are detected but
not supported; decompress them first with ``zstd -d``.

Requirements: numpy, scipy.
"""
import sys
import numpy as np
import scipy.sparse as sparse
import os
import gzip
import struct
import socket
from warnings import warn

class cell_ndarray(np.ndarray):
    '''numpy.ndarray subclass that can carry extra attributes such as ``header``.

    Args:
        input_array: anything accepted by numpy.asarray (not copied if possible).
        **attrs: attributes to attach to the new array.
    '''
    def __new__(cls, input_array, **attrs):
        obj = np.asarray(input_array).view(cls)
        obj.__dict__.update(attrs)
        return obj

    def __array_finalize__(self, obj):
        if obj is not None:
            self.__dict__.update(getattr(obj, "__dict__", {}))

# magic number are used in bin files to indicate data type (last 2 bytes; see sys/bin.h)
# dname describes data type. M_ types are fundamental types.
# MC_ types has been superseded by MCC_ANY
# Derived types (M_LOC64=0x036402, M_MAP64=0x016402, ...) have extra bits above
# the last 2 bytes and are read as the underlying dense type.
magic2dname={
    0x6400: 'M_CSP64',
    0x6401: 'M_DSP64',
    0x6402: 'M_DBL',
    0x6403: 'M_INT64',
    0x6404: 'M_CMP',
    0x6405: 'M_INT32',
    0x6406: 'M_CSP32', #32 is for integer index
    0x6407: 'M_DSP32',
    0x6408: 'M_FLT',
    0x6409: 'M_ZMP',
    0x640A: 'M_UINT8',
    0x640B: 'M_INT16',
    #0x6410: 'MC_CSP', #obsolete
    #0x6411: 'MC_SP',#obsolete
    #0x6412: 'MC_DBL',#obsolete
    #0x6413: 'MC_INT64',#obsolete
    #0x6414: 'MC_CMP',#obsolete
    #0x6415: 'MC_INT32',#obsolete
    0x6421: 'MCC_ANY',
    0x6430: 'M_SSP64',
    0x6431: 'M_SSP32',
    0x6432: 'M_ZSP64',
    0x6433: 'M_ZSP32',
}
M_0=0x0        #undefined. Used for empty elements of a cell array whose elements are all empty.
M_EOD=0x64FF   #end of data marker (used in socket streams)
M_SKIP=0x6600  #padding before the magic number
M_COMMENT=0x6500 #header (keywords) block
dname2type={#numpy data type of the values
    'M_ZMP': np.complex64,
    'M_CMP': np.complex128,
    'M_DBL': np.double,
    'M_FLT': np.float32,

    'M_UINT8': np.uint8,
    'M_INT16': np.int16,
    'M_INT32': np.int32,
    'M_INT64': np.int64,

    'M_CSP32': np.complex128,
    'M_CSP64': np.complex128,
    'M_DSP32': np.double,
    'M_DSP64': np.double,

    'M_ZSP32': np.complex64,
    'M_ZSP64': np.complex64,
    'M_SSP32': np.float32,
    'M_SSP64': np.float32
}

bitpix2magic={#FITS BITPIX to MAOS magic number for real values
    -32:0x6408,
    -64:0x6402,
    64:0x6403,
    32:0x6405,
    16:0x640B,
    8:0x640A,
    0:0
}
bitpix2magic_complex={ #for complex values (NAXIS=3 with NAXIS1=2)
    -32:0x6409,
    -64:0x6404
}
def readuint16(fp):
    '''Read a little endian uint16. Raises struct.error at end of file.'''
    return struct.unpack("<H", fp.read(2))[0]
def readuint32(fp):
    '''Read a little endian uint32. Raises struct.error at end of file.'''
    return struct.unpack("<I", fp.read(4))[0]
def readuint64(fp):
    '''Read a little endian uint64. Raises struct.error at end of file.'''
    return struct.unpack("<Q", fp.read(8))[0]
def readvec(fp, datatype, nxy):#enable read from socket and file
    '''Read nxy elements of numpy dtype datatype from fp.

    Works with regular files, gzip files, and sockets (file objects from
    open(fd)). Returns a writable 1-d array, which is shorter than nxy if the
    end of file is reached.
    '''
    #fp.fromefile is not good for gzip file.
    if fp.seekable() and fp.__class__!=gzip.GzipFile:
        return np.fromfile(fp, dtype=datatype, count=nxy, sep='')
    else:#for socket reading, or gzip
        #read into a preallocated array so that the result is writable (np.frombuffer is read only)
        out=np.empty(nxy, dtype=datatype)
        buf=memoryview(out.view(np.uint8))
        nread=0
        while nread<len(buf):
            nr=fp.readinto(buf[nread:])
            if not nr: #end of file
                break
            nread+=nr
        if nread<len(buf):
            out=out[0:nread//datatype.itemsize]
        return out

def convert_output(input, headers):
    '''Convert array to dictionary if necessary.

    If headers[0] is 'type=struct' and headers has one more entry than input,
    returns a dict mapping headers[i+1] to input[i]. Otherwise returns input unchanged.
    '''
    if len(headers)>0 and headers[0]=='type=struct' and len(headers)==input.size+1:
        out={}
        for ic in range(input.size):
            if len(headers)>ic+1:
                key=headers[ic+1]
            else:
                key='key_{}'.format(ic)
            while isinstance(key, list):
                key=key[0]
            out[key]=input[ic]
        return out
    else:
        return input

def readbin(file):
    '''Read a MAOS .bin or .fits file (optionally gzip compressed) or a socket.

    Args:
        file: file name (str or path-like; a leading ~ is expanded) or a
            connected socket.socket. If the file does not exist, the suffixes
            .bin, .bin.gz, .fits, .fits.gz (and .bin.zst, .fits.zst) are tried
            in turn, like MAOS does.

    Returns:
        The data: numpy array (cell_ndarray with a ``header`` attribute),
        scipy.sparse.csr_matrix, or numpy object array for cell arrays (see the
        module docstring for the shape conventions). A FITS file with several
        images returns an array of them.

        If the file does not exist, is zstd compressed, or cannot be parsed, a
        warning is issued and an empty array is returned (the file may not be
        ready yet, e.g., while a simulation is running).
    '''
    isfits=False
    isfile=False
    issock=False
    out=np.array(())
    err=0
    if isinstance(file, socket.socket):
        file=file.fileno()
        isfits=False
        isfile=False
        issock=True
    else:
        file=os.fspath(file)
        if file[0:1]=='~':
            file=os.path.expanduser(file)
        if not os.path.isfile(file): #file not found
            for suffix in ('.bin', '.bin.gz', '.fits', '.fits.gz', '.bin.zst', '.fits.zst'):
                if os.path.isfile(file+suffix):
                    file=file+suffix
                    break
            else:
                warn(f'readbin: file does not exist: {file}')
        if os.path.isfile(file):
            isfile=True
            if (file[-5:]=='.fits' or file[-8:] == '.fits.gz' or file[-9:] == '.fits.zst'):
                isfits=True
    if isfile or issock:
        gzipped=False
        if isfile:
            with open(file, 'rb') as fp:
                magic=fp.read(4)
            gzipped=magic[0:2]==b'\x1f\x8b'
            if magic==b'\x28\xb5\x2f\xfd':
                warn(f'readbin: {file} is zstd compressed, which is not supported. Decompress it first with zstd -d.')
                return out
        try:
            with (gzip.open(file, 'rb') if gzipped else open(file, 'rb', closefd=isfile)) as fp:
                (out, err)=readbin_auto(fp, isfits)
        except (OSError, EOFError, ValueError, KeyError, IndexError, struct.error) as error:#file may not be ready
            warn(f"readbin failed: {file}, {error}")
    if err<0:
        print('readbin: error happened', err)
    return out

def readbin_auto(fp, isfits, scan=0):
    '''Wraps readbin_do to handle top level scanning of fits file.

    Args:
        fp: opened binary file object.
        isfits: True for a FITS file.
        scan: if nonzero, read all the arrays until the end of file (always the
            case for FITS files, where each HDU is a separate array).

    Returns:
        (out, err): the data and the error code of readbin_do (0: success,
        1: end of data, -1: error). When scanning, err is 1 after all the data
        is read successfully. A single array is returned
        as is; several arrays are returned as an array of them (a numeric
        array if they have the same shape, otherwise an object array).
    '''
    err=0
    if scan or isfits:
        out=list()
        while err==0:
            [out_i, err]=readbin_do(fp, isfits)
            if err==0:
                out.append(out_i)
        if len(out)==1:
            out=out[0]
        else:
            try:
                out=np.array(out)
            except ValueError: #arrays of different shapes. Build the object array element by element.
                out2=np.empty(len(out), dtype=object)
                for i, out_i in enumerate(out):
                    out2[i]=out_i
                out=out2
    else:
        [out, err]=readbin_do(fp, isfits)

    return (out, err)

def readbin_do(fp, isfits):
    '''Read the next array (recursively for cell arrays) from fp.

    Args:
        fp: opened binary file object positioned at the start of a header.
        isfits: True for a FITS file.

    Returns:
        (out, err): out is the data (an empty array if there is none) and err
        is 0 on success, 1 at the end of data (M_EOD or end of file), or -1
        for an unknown data type.

    Raises:
        ValueError: if the magic number is invalid or the header is corrupted.
    '''
    err=0
    if isfits:
        [magic, nx, ny, header]=readfits_header(fp)
    else:
        [magic, nx, ny, header]=readbin_header(fp)

    if magic<0:
        dname='Unknown'
        err=-1
    elif magic==M_EOD:
        dname='EOD'
        err=1
    elif magic==M_0: #empty element in a cell array
        dname='M_0'
        nx=0
        ny=0
    else:
        try:
            dname=magic2dname[magic & 0xFFFF]
        except KeyError:
            raise ValueError("Invalid magic number {:#x} at {}".format(magic, fp.tell())) from None

    if dname[0:2]=='MC': #cell array
        out=np.zeros((ny, nx), dtype=object)

        for iy in range(0, ny):
            for ix in range(0, nx):
                out[iy,ix], err=readbin_do(fp, isfits)
                if err:
                    break
            if err:
                break
        if ny==1:
            out=out[0,]
    elif 'SP' in dname and nx>0 and ny>0: #sparse matrix
        #Layout: nz (uint64), column pointers (ny+1), row indices (nz), values (nz).
        #The index type is int64 for *SP64 and int32 for *SP32.
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
        if len(Jc)!=ny+1 or len(Ir)!=nz or len(P)!=nz:
            raise ValueError('readbin: sparse matrix is truncated')
        #MAOS CSC of nx*ny is python CSR of ny*nx
        out=sparse.csr_matrix((P, Ir, Jc), shape=(ny, nx))
    elif dname[0:2]=='M_' and nx>0 and ny>0: #dense matrix
        datatype=np.dtype(dname2type[dname])
        if isfits:
            datatype=datatype.newbyteorder('>')
        out=readvec(fp, datatype, nx*ny)
        if len(out) != nx*ny:
            warn('readbin: wrong length is read. Got {} expect {} at {} in {}'.format(len(out), nx*ny, fp.tell(), getattr(fp, 'name', '')))
        elif ny>1:
            out.shape=(ny, nx)

        if isfits:
            #data is padded to multiple of 2880 bytes
            byteread=datatype.itemsize*nx*ny
            byteleft=byteread%2880
            if byteleft:
                fp.read(2880-byteleft) #avoid seek, which doesn't work with gzip
            bscale=float(header.get('BSCALE', '1').split('/')[0])
            bzero=float(header.get('BZERO', '0').split('/')[0])
            if bscale!=1 or bzero!=0: #physical value
                out=out*bscale+bzero
                if datatype.kind in 'iu' and bscale==1 and bzero==int(bzero): #keep integer type, e.g., signed 8 bit
                    out=out.astype(np.int64 if bzero<0 or datatype.itemsize>4 else np.uint64)
    else:
        out=np.array(())
    if isinstance(out, np.ndarray):
        out=out.view(cell_ndarray)
    out.header=header #scipy sparse matrices also accept attributes
    return (out, err)


def readfits_header(fp):
    '''Read the header of the next FITS HDU (only IMAGE extensions are supported).

    Returns:
        (magic, nx, ny, header): MAOS magic number (M_EOD at the end of file),
        dimensions (NAXIS1, NAXIS2; for NAXIS=3, NAXIS1 must be 2 for complex
        numbers or 1), and a dict of the other keywords (values are the raw
        strings after '= ', including any comment).

    Raises:
        ValueError: for unsupported or corrupted headers.
    '''
    END=0
    page=0
    bitpix=0
    nx=0
    ny=0
    nc=1
    header={}
    while END==0:
        start=0
        if page==0: #first page
            res=fp.read(80)
            if len(res)!=80: #end of file
                END=2
                break
            elif res[0:6] !=b'SIMPLE' and res[0:16] != b"XTENSION= 'IMAGE":
                raise ValueError('Unknown data format in fits file at {}: {}'.format(fp.tell(), res))
            res=fp.read(80)
            bitpix=int(res[10:30])
            res=fp.read(80)
            naxis=int(res[10:30])
            if naxis>3:
                raise ValueError('readbin only support naxis<=3')
            if naxis>2:
                res=fp.read(80)
                nc=int(res[10:30])
                if nc!=2 and nc!=1:
                    raise ValueError('if naxis==3, axis[0] must be 2 for complex numbers')
            else:
                nc=1
            if naxis>0:
                res=fp.read(80)
                nx=int(res[10:30])
            if naxis>1:
                res=fp.read(80)
                ny=int(res[10:30])
            elif naxis==1: #a vector
                ny=1
            start=3+naxis
        for i in range(start, 36):
            res=fp.read(80)
            if len(res)!=80:
                raise ValueError('readbin: fits header is truncated')
            if not END:
                key=res[0:7].strip().decode('ascii')
                value=res[9:].strip().decode('ascii')
                if key=='END':
                    END=1 #do not break. let all data be read
                elif key=='COMMENT' or key == 'HISTORY':
                    header[key]=header.get(key,'')+value
                elif key=='PCOUNT':
                    value=int(value.split('/')[0])
                    if value != 0:
                        raise(ValueError('readbin only support PCOUNT=0, but is {}'.format(value)))
                elif key=='GCOUNT':
                    value=int(value.split('/')[0])
                    if value != 1:
                        raise(ValueError('readbin only support GCOUNT=1, but is {}'.format(value)))
                elif key!='EXTEND':
                    header[key]=value
        page=page+1
    if END == 1:
        try:
            if nc==2:#complex
                magic=bitpix2magic_complex[bitpix]
            else:#real
                magic=bitpix2magic[bitpix]
        except KeyError:
            raise ValueError('readbin: unsupported fits BITPIX={} (complex={})'.format(bitpix, nc==2)) from None
        return (magic, nx, ny, header)
    elif END == 2: #end of file
        return (M_EOD, 0, 0, header)
    else:
        return (-1, 0, 0, header)
def readbin_magic(fp):
    '''Read the next magic number, skipping M_SKIP padding. Returns M_EOD at end of file.'''
    try:
        magic=readuint32(fp)
        while magic==M_SKIP: #padding
            magic=readuint32(fp)
        return magic
    except struct.error: #end of file
        return M_EOD
def readbin_header(fp):
    '''Read the header of the next array in a .bin file.

    The layout (see write_bin_header in sys/bin.c) is zero or more keyword
    blocks [M_COMMENT, nlen (uint64), nlen bytes of text, nlen, M_COMMENT],
    an optional M_SKIP padding, then the magic number (uint32), nx and ny
    (uint64), all little endian. M_EOD (end of data) has no dimensions.

    Returns:
        (magic, nx, ny, header): header is a string, or a numpy array of
        strings if it has several lines.

    Raises:
        ValueError: if a keyword block is corrupted.
    '''
    header=''
    magic=readbin_magic(fp)
    while magic&0xFFFF==M_COMMENT:
        nlen=readuint64(fp)
        header=header+fp.read(nlen).decode('utf-8', errors='replace').strip('\x00')
        #header+=fp.read(nlen).strip() #.decode('utf-8')
        nlen2=readuint64(fp)
        magic2=readbin_magic(fp)
        if nlen!=nlen2 or magic!=magic2:
            raise ValueError('Header verification failed at {}'.format(fp.tell()))
        magic=readbin_magic(fp)
    if magic!=M_EOD:
        nx=readuint64(fp)
        ny=readuint64(fp)
    else:
        nx=0
        ny=0
    header=header.strip('\n')
    if header.find('\n')!=-1:
        header=np.array(header.split('\n'))
    return (magic, nx, ny, header)

if __name__ == '__main__':
    if len(sys.argv)>1:
        for fn in sys.argv[1:]:
            res=readbin(fn)
            print(fn, ' is ', res.shape)
    else:
        raise SystemExit("Usage: readbin.py file1 [file2 ...]\n   or in python: from readbin import readbin; res=readbin(file name)")
