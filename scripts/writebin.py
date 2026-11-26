#!/usr/bin/env python3
"""Pure Python writer for the MAOS .bin and .fits file formats.

This is the counterpart of readbin.py and writebin.m. It writes numpy arrays,
scipy sparse matrices, and (nested) cell arrays in the format that MAOS
(sys/bin.c, math/cell.c, math/spbin.c) and readbin.py/readbin.m read, without
needing the compiled aolib.so (see interface.py/libaos.py write() for the
ctypes based alternative).

Usage:
    from writebin import writebin
    writebin(opd, 'opd.bin', 'dx=1/64\\nox=-0.5\\noy=-0.5')
    writebin([np.eye(3), sparse.eye(3)], 'cell.bin.gz')
    python3 writebin.py in.bin out.fits  #convert a file read by readbin

Conventions (the same as readbin.py, so data read by readbin can be written back):
    - A numpy array of shape (ny, nx) is a MAOS array of nx*ny (column major);
      a 1-d array of shape (nx,) is nx*1. A 0-d array or scalar is 1*1.
    - A cell array is a numpy object array, list, or tuple (one that numpy
      cannot convert to a numeric array, or any object array). Elements may be
      None, which is written as an empty array like MAOS does.
    - A scipy sparse matrix of shape (ny, nx) is MAOS CSC of nx*ny. It is
      written with 64 bit indices.
    - The header (keywords) is taken from the ``header`` argument, or from the
      ``header`` attribute of the data (as set by readbin) if not given.

File name: if it does not end with .bin, .fits, .bin.gz, .fits.gz, .bin.zst,
or .fits.zst, .bin is appended, like MAOS does. .gz files are compressed with
gzip. .zst files require the zstandard module. .fits files support dense
arrays and cells of them only (sparse matrices are not supported by FITS).

Requirements: numpy, scipy (for sparse matrices).
"""
import sys
import os
import gzip
import socket
import struct
import numpy as np
try:
    import scipy.sparse as sparse
except ImportError: #sparse matrices are then not supported
    sparse=None

#magic numbers (see sys/bin.h and readbin.py)
M_CSP64=0x6400
M_DSP64=0x6401
M_DBL=0x6402
M_INT64=0x6403
M_CMP=0x6404
M_INT32=0x6405
M_FLT=0x6408
M_ZMP=0x6409
M_UINT8=0x640A
M_INT16=0x640B
MCC_ANY=0x6421
M_SSP64=0x6430
M_ZSP64=0x6432
M_COMMENT=0x6500
M_SKIP=0x6600

dtype2magic={#dense numpy data types written as is
    np.dtype(np.float64): M_DBL,
    np.dtype(np.float32): M_FLT,
    np.dtype(np.complex128): M_CMP,
    np.dtype(np.complex64): M_ZMP,
    np.dtype(np.int64): M_INT64,
    np.dtype(np.int32): M_INT32,
    np.dtype(np.int16): M_INT16,
    np.dtype(np.uint8): M_UINT8,
}
spdtype2magic={#sparse data types
    np.dtype(np.float64): M_DSP64,
    np.dtype(np.float32): M_SSP64,
    np.dtype(np.complex128): M_CSP64,
    np.dtype(np.complex64): M_ZSP64,
}
magic2bitpix={#FITS BITPIX. Complex numbers have an extra first axis of 2.
    M_FLT: -32, M_ZMP: -32, M_DBL: -64, M_CMP: -64,
    M_INT64: 64, M_INT32: 32, M_INT16: 16, M_UINT8: 8,
}
FITS_BS=2880 #FITS block size

def convert_dense(arr):
    '''Convert a numeric numpy array to a data type supported by MAOS.

    bool becomes uint8, int8 becomes int16, uint16 and uint32 the next larger
    signed integer, uint64 becomes int64 (values must fit), float16 becomes float32,
    and long double (complex) becomes double (complex).

    Raises:
        TypeError: for non-numeric data types.
        OverflowError: for uint64 values that do not fit in int64.
    '''
    dt=arr.dtype
    if dt in dtype2magic:
        return arr
    if dt.kind=='b':
        return arr.astype(np.uint8)
    if dt.kind=='u':
        if dt.itemsize==8:
            if arr.size and arr.max()>np.iinfo(np.int64).max:
                raise OverflowError('writebin: uint64 value does not fit in int64')
            return arr.astype(np.int64)
        return arr.astype({2:np.int32, 4:np.int64}[dt.itemsize])
    if dt.kind=='i': #int8 and other integer sizes, if any
        return arr.astype(np.int16 if dt.itemsize==1 else np.int64)
    if dt.kind=='f':
        return arr.astype(np.float32 if dt.itemsize<4 else np.float64)
    if dt.kind=='c':
        return arr.astype(np.complex128)
    raise TypeError(f'writebin: unsupported data type {dt}')

def issparse(data):
    '''Is data a scipy sparse matrix or array.'''
    return sparse is not None and sparse.issparse(data)

def to_array(data):
    '''Convert data to a numpy array (numeric, or object for cell arrays).

    None becomes an empty cell. A list or tuple becomes a 1-d object array
    (cell) if any element is a numpy array, sparse matrix, or None, or if numpy
    cannot convert it to a numeric array of at most 2 dimensions. Otherwise,
    e.g. [1, 2] or [[1, 2], [3, 4]], it becomes a numeric array.
    '''
    if data is None:
        return np.empty((0,), dtype=object)
    if isinstance(data, np.ndarray):
        return data
    if isinstance(data, (list, tuple)):
        if not any(x is None or issparse(x) or isinstance(x, np.ndarray) for x in data):
            try:
                arr=np.asarray(data)
                if arr.dtype!=object and arr.ndim<=2:
                    return arr
            except ValueError: #ragged
                pass
        out=np.empty(len(data), dtype=object)
        for i, x in enumerate(data):
            out[i]=x
        return out
    if isinstance(data, (str, bytes, dict)):
        raise TypeError(f'writebin: unsupported data type {type(data)}')
    return np.asarray(data)

def dims(arr):
    '''Return MAOS (nx, ny) of a numpy array: shape (ny, nx), (nx,), or ().

    Raises:
        ValueError: for arrays with more than 2 dimensions.
    '''
    if arr.ndim==0:
        return (1, 1)
    elif arr.ndim==1:
        return (arr.shape[0], 1)
    elif arr.ndim==2:
        return (arr.shape[1], arr.shape[0])
    else:
        raise ValueError(f'writebin: array of {arr.ndim} dimensions is not supported (at most 2)')

def header2str(header):
    '''Convert a header to a string, or None if empty.

    A numpy array or list of strings (multi-line header from readbin) is joined
    by newlines. A dict (FITS keywords from readbin) becomes key=value lines,
    with COMMENT and HISTORY values as plain text lines.
    '''
    if header is None:
        return None
    if isinstance(header, dict):
        lines=[]
        for key, val in header.items():
            if key in ('COMMENT', 'HISTORY'):
                lines.append(str(val))
            else:
                lines.append(f'{key}={val}')
        header='\n'.join(lines)
    elif isinstance(header, (np.ndarray, list, tuple)):
        header='\n'.join(str(h) for h in np.ravel(np.asarray(header, dtype=object)))
    elif isinstance(header, bytes):
        header=header.decode('utf-8')
    else:
        header=str(header)
    return header if header else None

def magic_of(data):
    '''Magic number of the element type of data (MCC_ANY for cells), used for
    the empty elements of a cell array, like MAOS does. 0 if unknown.'''
    if data is None:
        return 0
    if issparse(data):
        return spdtype2magic.get(np.dtype(data.dtype), M_DSP64)
    try:
        arr=to_array(data)
    except TypeError:
        return 0
    if arr.dtype==object:
        return MCC_ANY
    try:
        return dtype2magic[convert_dense(arr).dtype]
    except (TypeError, OverflowError):
        return 0

def fits_cards(text):
    '''Convert a header string to FITS cards like write_fits_header in sys/bin.c.

    Each key=value entry, separated by newline or ;, becomes a card KEY = VALUE
    (the key and the value before any / are converted to upper case). Other text
    becomes COMMENT cards. Text longer than 70 characters is continued in
    CONTINUE cards.

    Returns:
        list of 80 character strings.
    '''
    cards=[]
    if not text:
        return cards
    s=text.lstrip()
    cont=False
    upper=False
    while s:
        card=[' ']*80
        nl=len(s)
        for sep in ('\n', ';'):
            i=s.find(sep)
            if i!=-1 and i<nl:
                nl=i
        eq=-1 if cont else s.find('=', 0, nl)
        if eq!=-1:
            for ic in range(min(8, eq)):
                card[ic]=s[ic].upper()
            card[8]='='
            start=eq+1
            upper=True
        else:
            card[0:8]=list('CONTINUE' if cont else 'COMMENT ')
            start=0
            if not cont:
                upper=False
        length=nl-start
        truncated=False
        if length>70: #each card can contain maximum 70 characters
            length=70
            while length>1 and s[start+length-1].isspace(): #do not split at a space
                length-=1
            truncated=True
        for ic in range(length):
            c=s[start+ic]
            if c=='/':
                upper=False
            card[10+ic]=c.upper() if upper else c
        cards.append(''.join(card))
        if truncated:
            s=s[start+length:]
            cont=True
        else:
            s=s[nl+1:].lstrip() if nl<len(s) else ''
            cont=False
    return cards

class Writer:
    '''Writes arrays to an opened binary file object in .bin or .fits format.'''
    def __init__(self, fp, isfits):
        self.fp=fp
        self.isfits=isfits
        self.nhdu=0 #number of FITS HDUs written
        self.gstr=None #FITS: header of the enclosing cell, used for elements without header

    def write_header(self, magic, nx, ny, text):
        '''Write the header of an array (see write_header in sys/bin.c).'''
        if self.isfits:
            if magic==MCC_ANY: #cells are flattened. Keep the header for the elements.
                if text:
                    self.gstr=text
            else:
                self.write_fits_header(magic, nx, ny, text if text else self.gstr)
        else:
            buf=b''
            if text:
                #M_COMMENT, length, text padded with NUL to a multiple of 8 bytes (at least one NUL), length, M_COMMENT
                str8=text.encode('utf-8')
                nlen8=((len(str8)+8)>>3)<<3
                buf=struct.pack('<IQ', M_COMMENT, nlen8)+str8.ljust(nlen8, b'\0')+struct.pack('<QI', nlen8, M_COMMENT)
            self.fp.write(buf+struct.pack('<IIQQ', M_SKIP, magic, nx, ny))

    def write_fits_header(self, magic, nx, ny, text):
        '''Write a FITS HDU header (see write_fits_header in sys/bin.c).'''
        if magic not in magic2bitpix:
            raise TypeError(f'writebin: data type {magic:#x} (e.g. sparse) is not supported in fits')
        naxes=[nx, ny] if nx and ny else []
        if naxes and magic in (M_CMP, M_ZMP):
            naxes=[2]+naxes
        if self.nhdu==0:
            cards=['%-8s= %20s'%('SIMPLE', 'T')]
        else:
            cards=['%-8s= %s'%('XTENSION', "'IMAGE   '")]
        cards.append('%-8s= %20d'%('BITPIX', magic2bitpix[magic]))
        cards.append('%-8s= %20d'%('NAXIS', len(naxes)))
        for i, n in enumerate(naxes):
            cards.append('%-5s%-3d= %20d'%('NAXIS', i+1, n))
        if self.nhdu==0:
            cards.append('%-8s= %20s'%('EXTEND', 'T'))
        else:
            cards.append('%-8s= %20s'%('PCOUNT', '0'))
            cards.append('%-8s= %20s'%('GCOUNT', '1'))
        cards+=fits_cards(text)
        cards.append('END')
        cards+=['']*((-len(cards))%36) #pad to full pages of 36 cards
        self.fp.write(''.join(c.ljust(80) for c in cards).encode('ascii', errors='replace'))
        self.nhdu+=1

    def write_dense(self, arr, text):
        '''Write a numeric numpy array.'''
        arr=convert_dense(arr)
        nx, ny=dims(arr)
        magic=dtype2magic[arr.dtype]
        if arr.size==0:
            nx=ny=0
        self.write_header(magic, nx, ny, text)
        if arr.size:
            if self.isfits: #big endian, padded with zeros to a full block
                buf=np.ascontiguousarray(arr, dtype=arr.dtype.newbyteorder('>')).tobytes()
                self.fp.write(buf+b'\0'*((-len(buf))%FITS_BS))
            else:
                self.fp.write(np.ascontiguousarray(arr, dtype=arr.dtype.newbyteorder('<')).tobytes())

    def write_sparse(self, data, text):
        '''Write a scipy sparse matrix (see X(spwritedata) in math/spbin.c).'''
        arr=data.tocsr() #python CSR of ny*nx is MAOS CSC of nx*ny
        if arr.dtype not in spdtype2magic:
            arr=arr.astype(np.complex128 if arr.dtype.kind=='c' else np.float64)
        magic=spdtype2magic[np.dtype(arr.dtype)]
        if self.isfits:
            raise TypeError('writebin: sparse matrix is not supported in fits')
        nz=int(arr.indptr[-1]) if arr.shape[0] else 0
        if nz==0: #MAOS writes empty sparse matrix without dimension
            self.write_header(magic, 0, 0, None)
            return
        arr=arr.copy() #do not modify the input when sorting
        arr.sort_indices()
        ny, nx=arr.shape
        self.write_header(magic, nx, ny, text)
        self.fp.write(struct.pack('<Q', nz))
        self.fp.write(arr.indptr.astype('<i8').tobytes())
        self.fp.write(arr.indices[:nz].astype('<i8').tobytes())
        self.fp.write(np.ascontiguousarray(arr.data[:nz], dtype=arr.dtype.newbyteorder('<')).tobytes())

    def write(self, data, header=None):
        '''Write data (recursively for cells) with header (see writebin()).'''
        if header is None:
            header=getattr(data, 'header', None)
        if issparse(data):
            self.write_sparse(data, header2str(header))
            return
        arr=to_array(data)
        if arr.dtype!=object:
            self.write_dense(arr, header2str(header))
            return
        #cell array (see writedata in math/cell.c)
        nx, ny=dims(arr)
        nel=arr.size
        if nel==0:
            nx=ny=0
        hcell=header
        helem=None
        if isinstance(header, (list, tuple)) and nel>0 and len(header) in (nel, nel+1):
            #one header per element, optionally followed by that of the cell (as in writebin.m)
            helem=header[0:nel]
            hcell=header[nel] if len(header)==nel+1 else None
        saved_gstr=self.gstr
        self.write_header(MCC_ANY, nx, ny, header2str(hcell))
        elems=arr.ravel(order='C') #(ny, nx) row major is column major nx*ny
        magic0=0 #type of empty elements: that of the first non-empty element
        for x in elems:
            if x is not None:
                magic0=magic_of(x)
                break
        for i, x in enumerate(elems):
            hi=helem[i] if helem is not None else None
            if x is None:
                if self.isfits: #an image without data
                    self.write_fits_header(magic0 if magic0 in magic2bitpix else M_DBL, 0, 0, header2str(hi) or self.gstr)
                else:
                    self.write_header(magic0, 0, 0, header2str(hi))
            else:
                self.write(x, hi)
        self.gstr=saved_gstr

def outname(fn):
    '''Expand ~ and append .bin if fn does not have a recognized suffix, like MAOS does.'''
    fn=os.path.expanduser(os.fspath(fn))
    if not fn.endswith(('.bin', '.fits', '.bin.gz', '.fits.gz', '.bin.zst', '.fits.zst')):
        fn+='.bin'
    return fn

def writebin(data, file, header=None):
    '''Write data to a MAOS .bin or .fits file, or a socket.

    Args:
        data: numpy array (at most 2-d), scalar, scipy sparse matrix, None, or
            a cell array: numpy object array, list or tuple of these (nested
            cells allowed). See the module docstring for the shape convention.
        file: file name (str or path-like; a leading ~ is expanded; .bin is
            appended if it has no recognized suffix; .gz/.zst compresses the
            file), or a connected socket.socket to stream the data in .bin
            format.
        header: header text, normally key=value entries separated by newlines
            (MAOS keywords, e.g., 'dx=0.1\\nox=-5'), or a list/array of lines,
            or a dict of FITS keywords. Defaults to the ``header`` attribute
            of data (as set by readbin), if any. For a cell array it may also
            be a list (or tuple, not a numpy array) with one header per
            element, optionally followed by the header of the cell itself
            (numel+1 entries). Element headers otherwise come from the
            ``header`` attribute of each element.

    Returns:
        The file name written (None for a socket).

    Raises:
        TypeError: for unsupported data types (e.g., str, dict, or sparse
            matrices in FITS files).
        ValueError: for arrays with more than 2 dimensions.
        ImportError: for .zst files if the zstandard module is not available.
        OSError: if the file cannot be written.
    '''
    if isinstance(file, socket.socket):
        with file.makefile('wb') as fp:
            Writer(fp, False).write(data, header)
            fp.flush()
        return None
    fn=outname(file)
    isfits=fn.endswith(('.fits', '.fits.gz', '.fits.zst'))
    #Write to a temporary file in the same folder and rename it at the end, so
    #that a reader never sees a partial file. This also replaces a symbolic
    #link instead of writing through it, like MAOS does.
    fntmp=os.path.join(os.path.dirname(fn) or '.', '.'+os.path.basename(fn)+'.%d.tmp'%os.getpid())
    try:
        if fn.endswith('.gz'):
            fp=gzip.open(fntmp, 'wb')
        elif fn.endswith('.zst'):
            try:
                import zstandard
            except ImportError:
                raise ImportError('writebin: the zstandard module is required to write '+fn) from None
            fp=zstandard.open(fntmp, 'wb')
        else:
            fp=open(fntmp, 'wb')
        with fp:
            Writer(fp, isfits).write(data, header)
        os.replace(fntmp, fn)
    finally:
        if os.path.exists(fntmp):
            os.remove(fntmp)
    return fn

if __name__ == '__main__':
    if len(sys.argv)==3:
        from readbin import readbin
        res=readbin(sys.argv[1])
        print('wrote', writebin(res, sys.argv[2]))
    else:
        raise SystemExit("Usage: writebin.py input output  #convert a file read by readbin\n"
                         "   or in python: from writebin import writebin; writebin(data, file name, header)")
