#!/Usr/bin/env python

#Use ctypes to interface with C libraries
#POINTER is the class type of pointer
#pointer() acts on actual array, while POINTER() works on class type.
#pointer(cell(arr)) #creates cell Structure for np.array type arr and makes a pointer
#pcell=POINTER(cell) ; pcell() #Creates a class for cell Structure and then makes an object pointer with no content.

#c_int is a ctypes type
#a=c_int(42) #creates a ctypes int object
#p=pointer(a) creates C compatible pointers to object a
#p.contents retreates the contents of pointer p
#addressof(p) retreates the address of p
#use pd=cast(p, POINTER(c_double)) to convert pointer p to c_double pointer 
#use pd=cast(address, POINTER(c_double)) to convert address in c_double pointer
#cast(addressof(a), POINTER(c_double)).contents #reinterpreted int a as double

#pointer() creates a real pointer to a ctype object (Sturecture)
#byref() is a simplified version of pointer(). It cannot be used as byref(byref())
#can use byref(pointer())
#Set restype to return correct value
#Set argtypes for type checking for input into C code
#be careful regarding memory management.


#TODO: investigate whether CFFI is a better solution.
import os
#import sys
#from pdb import set_trace as keyboard
from ctypes import *
#import json
import numpy as np
import scipy.sparse as sp
#from copy import deepcopy
#from warnings import warn
aolib_so=os.environ.get('MAOS_AOLIB', 'aolib.so')
try:
    lib=cdll.LoadLibrary(aolib_so)
except:
    raise Exception('Load aolib.so failed from '+aolib_so)
from readbin import set_header, convert_output

if 'c_double_complex' in globals():
    has_complex=1
else:
    has_complex=0
    class c_double_complex_struct(Structure):
        _fields_ = [
            ("real", c_double),
            ("imag", c_double)
        ]
    # Assign it to a usable name
    c_double_complex = c_double_complex_struct
    class c_float_complex_struct(Structure):
        _fields_ = [
            ("real", c_float),
            ("imag", c_float)
        ]
    # Assign it to a usable name
    c_float_complex = c_float_complex_struct
    
def simplify(arr, do_stack=1, do_squeeze=0):
    '''convert object array a[i,j][k,n]... to simple ndarray a[i,j,k,n,...]
        shape and singleton dimensions are preserved.
    '''
    if type(arr) is list:
        if len(arr)==1:
            arr=arr[0]
        else:
            arr=np.array(arr)
    if type(arr) is np.ndarray and do_squeeze:
        arr=np.squeeze(arr)
        if arr.size==1:
            return arr.item(0) #this works for 0d array
    if type(arr) is np.ndarray and arr.dtype.name=='object':
        flags=np.zeros(arr.shape,dtype=bool)
        can_stack=0
        for ind,iarr in np.ndenumerate(arr):
            arr[ind]=simplify(iarr,do_stack,do_squeeze)
            if arr[ind].size>0:
                flags[ind]=True
                if arr[ind].dtype.name!='object':
                    can_stack=1
                    shape=arr.shape+arr[ind].shape
        if do_stack and can_stack:
            try:
                arr=np.stack(arr[flags])
                arr.shape=shape
                if do_squeeze:
                    arr=np.queeze(arr)
            except Exception as err:
                #print('stack failed', err)
                pass
    return arr

#obtain ctypes type information from MAOS id.
#The value is (type, complex?, kind(0:dense, 1: sparse, 2: loc, 10: cell))
id2ctype={
    0x6400: (c_double_complex,0,1) if has_complex else (c_double,1,1), #M_CSP64
    0x6401: (c_double,0,1), #M_SP64
    0x6402: (c_double,0,0), #'M_DBL'
    0x6403: (c_long,  0,0), #'M_INT64'
    0x6404: (c_double_complex, 0,0) if has_complex else (c_double,1,0), #'M_CMP'
    0x6405: (c_int,   0,0), #'M_INT32',),
    0x6406: (c_float_complex, 0, 1) if has_complex else (c_float, 1,1),#'M_CSP32,
    0x6407: (c_float, 0,1), #'M_SP32',),
    0x6408: (c_float, 0,0), #'M_FLT',),
    0x6409: (c_float_complex, 0, 0) if has_complex else (c_float, 1,0), #'M_ZMP'
    0x640A: (c_char,  0,0), #'M_INT8',),
    0x640B: (c_short, 0,0), # 'M_INT16',),
    0x6421: (c_void_p,0,10),#MC_ANY
    0x36402: (c_double,0,2),#M_LOC64
}
id2dtype={
    0x6400: np.complex128,
    0x6401: np.float64,
    0x6402: np.float64,
    0x6403: np.int64,
    0x6404: np.complex128,
    0x6405: np.int32,
    0x6406: np.complex64,
    0x6407: np.float32,
    0x6408: np.float32,
    0x6409: np.complex64,
    0x640A: np.int8,
    0x640B: np.int16,
}
dtype2id={#Conversion from numpy dtype to maos id
    np.double:  0x6402,
    np.int64:   0x6403,
    np.complex128:0x6404,
    np.int32:   0x6405,
    np.float32: 0x6408,
    np.complex64:0x6409,
    np.int8:    0x640A,
    np.int16:   0x640B,
    np.object_: 0x6421
}
spdtype2id={#Conversion from sparse type to maos id
    np.float32:   0x6407,
    np.float64:   0x6401,
    np.complex64: 0x6406,
    np.complex128:0x6400,
}
def py2cell(arr, tid=0):
    '''convert numpy.array to C array adaptively (based on its type)'''
    if arr is None:
        return None
    if type(arr) is list:
        arr=np.asarray(arr)
    if sp.issparse(arr):
        return csr(arr, tid)
    else:
        return cell(arr, tid)

def py2cellref(arr, tid=0):
    '''convert numpy array to C array pointer adaptively or output'''
    if arr is None:
        return None
    elif type(arr) is list:
        arr = np.asarray(arr)
    if sp.issparse(arr):
        return byref(csr(arr,tid))
    elif type(arr) is np.ndarray:
        if arr.size==0:
            return None #turn empty ndarray to Null pointer. do not use 0
        else:
            return byref(cell(arr,tid))
    else:
        return byref(arr)

def arr2object(arr): 
    '''convert numpy number arrays with more than 3 dimensions to numpy object arrays'''
    if type(arr) is not np.ndarray:
        print('Unsupported data type')
        return arr
    if arr.ndim<=2:
        return arr
    ndim=2-(arr.ndim%2)
    arr2=np.empty((arr.shape[0:ndim]),dtype=object)
    if ndim==2:
        for i in range(arr2.shape[0]):
            for j in range(arr2.shape[1]):
                arr2[i,j]=arr2object(arr[i,j])
    else:
        for i in range(arr2.shape[0]):
            arr2[i]=arr2object(arr[i])
    return arr2

class wrap_pointer:
    '''Reference counting for c array'''
    def __init__(self, pointer):
        self.pointer=pointer
    def __del__(self):
        lib.cellfree_do(cast(self.pointer, c_void_p))
class cell_ndarray(np.ndarray):
    '''Subclass to manage memory and extra attributes'''
    def __new__(cls, ctypes_pointer, pointer=None):
        #print(f'cell_ndarray __new__ {ctypes_pointer}')
        ctypes_array=cast(ctypes_pointer, POINTER(cell)).contents
        if pointer is None and not hasattr(ctypes_array, 'python'):
            pointer=wrap_pointer(ctypes_pointer)
        if ctypes_array.header:
            header=ctypes_array.header.decode('ascii')
        else:
            header=''
        try:
            (tt, iscomplex, kind)=id2ctype.get(ctypes_array.id)
        except:
            print("id2ctype: unknown type", id)
            res=np.array([])
        if kind==0 or kind==2: #dense matrix
            parr=cast(ctypes_array.p, POINTER(tt))
            if iscomplex:
                res=np.ctypeslib.as_array(parr, shape=(*ctypes_array.shape(0),2))
                if tt is c_double:
                    res=np.squeeze(res.view(np.complex128), axis=-1)
                elif tt is c_float:
                    res=np.squeeze(res.view(np.complex64), axis=-1)
                else:
                    raise(Exception('Please implement'))
            else:
                res=np.ctypeslib.as_array(parr, shape=ctypes_array.shape(0))
            if kind==2: #LOC
                res.dx=ctypes_array.dx
                res.dy=ctypes_array.dy
                res.ht=ctypes_array.ht
                res.iac=ctypes_array.iac

        elif kind==1: #sparse matrix does not support subclassing like numpy
            return cell_csr_array(ctypes_pointer, pointer)
        elif kind==10: #cell array
            res=np.empty(ctypes_array.shape(1), dtype=object)
            parr=cast(ctypes_array.p, POINTER(c_void_p))
            for iy in range(ctypes_array.ny):
                for ix in range(ctypes_array.nx):
                    address=parr[ix+ctypes_array.nx*iy]
                    if address is not None:
                        pp=cast(int(address), POINTER(cell))
                        res[iy, ix]=cell_ndarray(pp, pointer) #recursive. do not keep pointer in child
                    else:
                        res[iy, ix]=np.array([])
            if ctypes_array.ny==1:
                res=res[0,]

        obj=res.view(cls)
        obj.pointer=pointer #pointer is reference counted
        obj.header=header
        if len(header)>0:
            for pair in header.replace('\n',';').replace(';;',';').split(';'):
                if pair.find('=')>-1:
                    ires=pair.split('=')
                    if len(ires)==2:
                        try:
                            setattr(obj, ires[0], float(ires[1]))
                        except:
                            setattr(obj, ires[0], ires[1])
        return obj
    def __array_finalize__(self, obj):
        if obj is None:
            return
        if hasattr(obj, 'pointer'):
            self.pointer=obj.pointer

class cell_csr_array(sp.csr_array):
    def __init__(self, ctypes_pointer, pointer=None):
        #print(f'cell_csr_array __new__ {ctypes_pointer}')
        ctypes_array=cast(ctypes_pointer, POINTER(csr)).contents
        if pointer is None and not hasattr(ctypes_array, 'python'):
            pointer=wrap_pointer(ctypes_pointer)
        try:
            (tt, iscomplex, kind)=id2ctype.get(ctypes_array.id)
        except:
            print("id2ctype: unknown type", id);
            return None
        if kind==0 or kind==2:
            return cell_ndarray(ctypes_array, pointer)
        elif kind==1:
            if ctypes_array.nzmax>0:
                xp=np.ctypeslib.as_array(cast(ctypes_array.x, POINTER(tt)), shape=(ctypes_array.nzmax,))
                ip=np.ctypeslib.as_array(cast(ctypes_array.i, POINTER(c_long)), shape=(ctypes_array.nzmax,))
                pp=np.ctypeslib.as_array(cast(ctypes_array.p, POINTER(c_long)), shape=(ctypes_array.ny+1,))
                super().__init__((xp, ip, pp), shape=(ctypes_array.ny, ctypes_array.nx), copy=False)
                self.pointer=pointer 
            else:
                super().__init__((ctypes_array.nx,ctypes_array.ny))
    def __array_finalize__(self, obj):
        if obj is None:
            return
        if hasattr(obj, 'pointer'):
            self.pointer=obj.pointer
            
def pt2py(pointer):
    '''convert C array pointer to numpy array. References C memory'''
    if bool(pointer):
        try:
            (tt, iscomplex, kind)=id2ctype.get(pointer.contents.id)
        except:
            print("id2ctype: unknown type", pointer.contents.id)
            return np.array([])
        if kind==1:#sparse
            return cell_csr_array(pointer)
        else:
            return cell_ndarray(pointer)
    else:
        return np.array([])

class cell(Structure):
    '''To interface numpy.array with C cell '''
    _fields_ = [ #fields compatible with C type of cell and mat
        ('id', c_uint32),
        ('p',  c_void_p),
        ('nx', c_long),
        ('ny', c_long),
        ('header', c_char_p),
        ('dummy_fp', c_void_p),
        ('dummy_fft', c_void_p),
        ('dummy_mem', c_void_p),
        ('dummy_async', c_void_p),
        ('dummy_deinit', c_void_p),
        ('dummy_make_keywords', c_void_p),
        ]
 
    def __init__(self, arr=None, tid=0):#convert from numpy to C. Memory is borrowed
        #attributes set within __init__ are per object
        #print(f'__init__ is called for cell {addressof(self)}')
        if type(arr) is list:
            arr=np.asarray(arr)
        if arr.ndim>2:
            arr=arr2object(arr)
        if arr is not None:
            tmpid=dtype2id.get(arr.dtype.type)
            if tmpid is None:
                print("init: Unknown data" +str( arr.dtype.type))
                return None
            if tid!=0 and tmpid != tid and tid!=0x6421 : #and tmpid !=0x6421:
                dtp = next(key for key, value in dtype2id.items() if value == tid) 
                try: #convert data
                    arr=arr.astype(dtp)
                    tmpid=tid
                except:
                    raise(Exception('data mismatch want {}, got {}'.format(tmpid, tid)))
            self.id=tmpid

            if arr.ndim>0:
                self.nx=arr.shape[-1]
            if arr.ndim>1:
                self.ny=arr.shape[-2]
            else:
                if self.nx>0:
                    self.ny=1
                else:
                    self.ny=0
            if self.nx==0:
                self.p=0
            elif arr.dtype.kind != 'O':
                if arr.flags['C']==False:
                    arr=arr.copy(order='C')
                self.p=arr.ctypes.data_as(c_void_p)
            else:
                self.qarr=np.zeros(self.shape(1), dtype=object) #stores objects
                self.parr=np.zeros(self.shape(1), dtype=c_void_p) #store pointer of objects
                for iy in range(self.ny):
                    for ix in range(self.nx):
                        if arr.ndim==1:
                            arri=arr[ix]
                        else:
                            arri=arr[iy,ix]
                        if arri is not None:
                            self.qarr[iy,ix]=py2cell(arri) #object array
                            self.parr[iy,ix]=addressof(self.qarr[iy,ix]) #object pointer array
                        else:
                            self.parr[iy,ix]=0
                self.p=self.parr.ctypes.data_as(c_void_p)
            self.arr=arr #keep a reference to arr as it can be a copy otherwise it may be destroyed
        else:
            self.id=0x6421
            self.p=None
            self.nx=0
            self.ny=0
        self.header=None
        self.dummy_fp=None
        self.dummy_fft=None
        self.dummy_mem=None
        self.dummy_async=None
        self.dummy_deinit=None
        self.dummy_make_keywords=None
        self.python=True
    def shape(self, twod):
        if self.ny > 1 or twod:
            return (self.ny, self.nx)
        else:
            return (self.nx,) #last , is necessary
            
class loc(Structure):
    '''To interface numpy.array with C lob '''
    _fields_ = [
        ('id',   c_uint32),
        ('locx', c_void_p),
        ('nloc', c_long),
        ('two',  c_long), 
        ('header',c_void_p),
        ('dummy_fp',c_void_p),
        ('dummy_fft',c_void_p),
        ('dummy_mem',c_void_p),
        ('dummy_async',c_void_p),
        ('dummy_deinit',c_void_p),
        ('dummy_make_keywords',c_void_p), #up to here. follow cell memory layout.
        ('locy', c_void_p),
        ('locstat_t', c_void_p),
        ('map',  c_void_p),
        ('dx',   c_double),
        ('dy',   c_double),
        ('ht',   c_double),
        ('iac',  c_double),
        ('dratio',c_double),
        ('npad', c_int),
        ]
    def __init__(self, arr=None): #convert from numpy to C. Memory is borrowed
        self.id= 0x36402 #M_LOC64
        if arr is not None:
            if len(arr.shape)!=2 or arr.shape[0] !=2 :
                raise(Exception('Array has to be of shape 2xn'))
            else:
                if arr.flags['C']==False:
                    arr=self.arr=arr.copy()
                if arr.dtype!=np.double:
                    arr=arr.astype(np.double)
                dlocx=np.abs(arr[0,1:]-arr[0,0:-1])
                dlocy=np.abs(arr[1,1:]-arr[1,0:-1])
                self.locx=arr[0,].ctypes.data_as(c_void_p)
                self.nloc=arr.shape[1]
                self.two=2;
                self.header=None
                self.dummy_fp=None
                self.dummy_fft=None
                self.dummy_mem=None
                self.dummy_async=None
                self.dummy_deinit=None
                self.dummy_make_keywords=None
                self.locy=arr[1,].ctypes.data_as(c_void_p)
                self.locstat_t=None
                self.map=None
                self.dx=min(dlocx[dlocx>0])
                self.dy=min(dlocy[dlocy>0])
                self.ht=0
                self.iac=0
                self.dratio=1
                self.npad=1
                self.python=True
                #print('loc: dx={0}, dy={1}'.format(self.dx, self.dy))
        #default initialization to zero

class csr(Structure):#CSR sparse matrix. We convert C CSC to Python CSR just like arrays as numpy is row order
    '''To interface numpy.array with C sparse array '''
    _fields_=[ #need to match C memory layout
        ('id', c_uint32),
        ('x', c_void_p),
        ('nx', c_long),
        ('ny', c_long),
        ('header',   c_char_p),
        ('dummy_fp', c_void_p),
        ('nzmax', c_long),
        ('p', c_void_p),
        ('i', c_void_p),
        ('nref', c_void_p),
    ]
    def __init__(self, arr=None, tid=0): #convert from numpy to C. Memory is borrowed
        if arr is None:
            self.id=spdtype2id.get(np.float64)
        elif sp.issparse(arr):
            if arr.dtype.type in (np.int32, np.int64):
                arr=arr.astype(np.float64)
            self.id=spdtype2id.get(arr.dtype.type)
            if  tid !=0 and tid !=0x6421 and self.id != tid:
                raise(Exception('data mismatch want {}, got {}'.format(self.id, tid)))
            if not arr.format in ('csr','csc'):
                print('Sparse format {} is not supported; will convert to csr'.format(arr.format))
                arr=arr.tocsr()
            if arr.format=='csr':
                self.ny, self.nx=arr.shape #python csr is maos csc transposed
            elif arr.format=='csc':
                self.nx, self.ny=arr.shape 
            #save subarrays to avoid being GC'ed.
            self.xp=arr.data
            self.ip=arr.indices.astype(np.long)
            self.pp=arr.indptr.astype(np.long) #p
            self.x=self.xp.ctypes.data_as(c_void_p) #data
            self.i=self.ip.ctypes.data_as(c_void_p) #row index
            self.p=self.pp.ctypes.data_as(c_void_p)
            self.nzmax=self.pp[-1]
            self.python=True
        else:
            raise(Exception('Invalid conversion to csr'))

def convert_fields(fields):
    '''convert a C type keyword to ctypes type'''
    val2type={
        '*':c_void_p,
        'double':c_double,
        'real':c_double,
        'float':c_float,
        'long':c_long,
        'int':c_int,
    }
    newfields=[]
    for key,val in fields.items():
        if val[-1]=='*':
            val=c_void_p
        else:
            val=val2type[val]
        newfields.append((key,val))
    return newfields

def make_class(name, fields):
    '''Dynaically create a ctypes class with field listed '''
    newfields=convert_fields(fields)
    class newclass(Structure):
        pass
        def as_array(self):#convert struct into dictionary
            out=dict()
            for ff in self._fields_:
                #convert C pointers to POINTER then to array
                field=eval("self.{}".format(ff[0]))
                if field is not None:
                    if ff[1] is c_void_p:
                        out[ff[0]]=cast(field,POINTER(cell)).contents.as_array()
                        #exec('out[\''+ff[0]+'\']=cast(self.'+ff[0]+',POINTER(cell)).contents.as_array()')
                    else:
                        out[ff[0]]=field
                else:
                    out[ff[0]]=None
            out['struct']=self
            return out

    newclass._fields_=newfields
    return newclass
