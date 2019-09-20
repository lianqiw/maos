#!/Usr/bin/env python

#Use ctypes to interface with C libraries
#POINTER is the class type of pointer
#c_int is a ctypes type
#a=c_int(42) #creates a ctypes int object
#p=pointer(a) creates pointers to object a
#p.contents retreates the contents of pointer p
#addressof(p) retreates the address of p
#use pd=cast(p, POINTER(c_double)) to convert pointer p to c_double pointer 
#use pd=cast(address, POINTER(c_double)) to convert address in c_double pointer
#cast(addressof(a), POINTER(c_double)).contents #reinterpreted int a as double

#Use byref() or pointer() to obtain pointer to an object
#Set restype to return correct value
#Set argtypes for type checking for input into C code
#be careful regarding memory management.


#TODO: investigate whether CFFI is a better solution.

import sys
from pdb import set_trace as keyboard
from ctypes import *
import json
import numpy as np
import scipy.sparse as sp
from warnings import warn

lib=cdll.LoadLibrary("/home/lianqiw/work/aos/comp/optim/bin/aolib.so")

id2ctype={
    #obtain type information from MAOS id.
    #The value is (type, is complex, kind(0:dense, 1: sparse, 2: loc, 10: cell))
    25600: (c_double,1,1), #M_CSP64
    25601: (c_double,0,1), #M_SP64
    25602: (c_double,0,0), #'M_DBL'
    25603: (c_long,  0,0), #'M_INT64'
    25604: (c_double,1,0), #'M_CMP'
    25605: (c_int,   0,0), #'M_INT32',),
    25606: (c_float, 1,1), #'M_CSP32',),
    25607: (c_float, 0,1), #'M_SP32',),
    25608: (c_float, 0,0), #'M_FLT',),
    25609: (c_float, 1,0), #'M_ZMP',),
    25610: (c_char,  0,0), #'M_INT8',),
    25611: (c_short, 0,0), # 'M_INT16',),
    25633: (c_void_p,0,10),#MC_ANY
    222210: (c_double,0,2),#M_LOC64
}
#convert C array pointer to numpy array. Freeing C memory
def pt2py(pointer):
    if bool(pointer):
        out=pointer.contents.as_array()
        pointer.contents.free()
        return out
    else:
        return None
#convert C vector to numpy array. Memory is copied.
def as_array(arr, id, shape):
    ''' convert C array arr to numpy based in id'''
    (tt, iscomplex, issparse)=id2ctype.get(id)
    if tt is None or not bool(arr) or shape[0]==0:
        return np.empty((0,))
    else:
        parr=cast(arr, POINTER(tt))
        if iscomplex:
            nparr=np.ctypeslib.as_array(parr, shape=(*shape,2))
            nparr2=nparr[...,0]+1j*nparr[...,1]
        else:
            nparr=np.ctypeslib.as_array(parr, shape=shape)
            nparr2=np.copy(nparr)
        return nparr2

#convert numpy array to any C array adaptively
def py2cell(arr):
    if type(arr) is list:
        arr=np.asarray(arr)
    if sp.isspmatrix_csc(arr):
        return csc(arr)
    else:
        return cell(arr)

#convert numpy array to any C array adaptively
def py2cellref(arr):
    if type(arr) is list:
        arr=np.asarray(arr)
    elif type(arr) is not np.ndarray:
        print('py2cellref does not take ', type(arr))
        return None
    if arr.size==0:
        return None #turn empty ndarray to Null pointer. do not use 0
    elif sp.isspmatrix_csc(arr):
        return byref(csc(arr))
    else:
        return byref(cell(arr))

class cell(Structure):
    _fields_ = [ #fields compatible with C type
        ('id', c_uint32),
        ('p',  c_void_p),
        ('nx', c_long),
        ('ny', c_long),
        ('header', c_char_p),
        ('mmap', c_void_p),
        ('nref', c_void_p),
        ('fft', c_void_p),
        ]
 
    def __init__(self, arr=None):#convert from numpy to C. Memory is borrowed
        dtype2id={#Conversion from numpy type to maos id
            np.double:25602,
            np.int64: 25603,
            np.object_:25633,
        }
        if type(arr) is list:
            arr=np.asarray(arr)

        if arr is not None:
            self.id=dtype2id.get(arr.dtype.type)
            if self.id is None:
                print("init: Unknown data" +str( arr.dtype.type))
                return None
            if arr.ndim>2:
                print("init: Only use 2 dimensions\n");
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
            elif arr.dtype.kind is not 'O':
                self.p=arr.ctypes.data_as(c_void_p)
            else:
                self.qarr=np.zeros(self.shape(1), dtype=object)
                self.parr=np.zeros(self.shape(1), dtype=c_void_p) #store pointers 
                for iy in range(self.ny):
                    for ix in range(self.nx):
                        if arr.ndim==1:
                            arri=arr[ix]
                        else:
                            arri=arr[iy,ix]
                        if arri is not None:
                            self.qarr[iy,ix]=py2cell(arri) #keep reference
                            self.parr[iy,ix]=addressof(self.qarr[iy,ix]) #pointer
                        else:
                            self.parr[iy,ix]=0
                self.p=self.parr.ctypes.data_as(c_void_p)
        else:
            self.id=25633
            self.p=None
            self.nx=0
            self.ny=0
        self.header=None
        self.mmap=None
        self.nref=None
        self.fft=None

    def shape(self, twod):
        if self.ny > 1 or twod:
            return (self.ny, self.nx)
        else:
            return (self.nx,) #last , is necessary

    def as_array(self): #convert form C to numpy. Memory is copied
        try:
            (tt, iscomplex, kind)=id2ctype.get(self.id)
        except:
            kind=-1

        if kind==0: #dense matrix
            if self.header:
                print(self.header)
            return as_array(self.p, self.id, self.shape(0))
        elif kind==1: #sparse matrix
            return cast(addressof(self), POINTER(csc)).contents.as_array()
        elif kind==2: #loc
            return cast(addressof(self), POINTER(loc)).contents.as_array()
        elif kind==10: #cell
            res=np.empty(self.shape(1), dtype=object)
            parr=cast(self.p, POINTER(c_void_p))
            for iy in range(self.ny):
                for ix in range(self.nx):
                    address=parr[ix+self.nx*iy]
                    if address is not None:
                        pp=cast(int(address), POINTER(cell))
                        res[iy, ix]=pp.contents.as_array() #recursive
                    else:
                        res[iy, ix]=np.empty(())
            if self.ny==1:
                res=res[0,]
            return res
        else:
            print('as_array: Unknown data, id='+ str(self.id))
            return np.empty((),dtype=object)
    def free(self):
        lib.cellfree_do(byref(self))

class loc(Structure):
    _fields_ = [
        ('id', c_uint32),
        ('locx',  c_void_p),
        ('locy',  c_void_p),
        ('nloc', c_long),
        ('dx', c_double),
        ('dy', c_double),
        ('ht', c_double),
        ('iac', c_double),
        ('locstat_t', c_void_p),
        ('map', c_void_p),
        ('npad', c_int),
        ('nref', c_void_p),
        ]
    def __init__(self, arr=None): #convert from numpy to C. Memory is borrowed
        self.id= 222210 #0x036402 #M_LOC64
        if arr is not None:
            if len(arr.shape)!=2 or arr.shape[0] !=2 :
                raise(Exception('Array has to of shape 2xn'))
            else:
                self.nloc=arr.shape[1]
                self.locx=arr[0,].ctypes.data_as(c_void_p)
                self.locy=arr[1,].ctypes.data_as(c_void_p)
                dlocx=arr[0,1:]-arr[0,0:-1]
                self.dx=min(dlocx[dlocx>0]);
                dlocy=arr[1,1:]-arr[1,0:-1]
                self.dy=min(dlocy[dlocy>0]);
                #print('loc: dx={0}, dy={1}'.format(self.dx, self.dy))
        else:
            self.nloc=0
            self.locx=None
            self.locy=None
            self.dx=0
            self.dy=0
        self.ht=0
        self.iac=0
        self.locstat_t=0
        self.map=0
        self.npad=0
        self.nref=0

    def as_array(self): #convert form C to numpy. Memory is copied
        if(self.locx):
            if self.id!=222210:
                raise(Exception('Wrong type'))
            else:
                arr=np.empty((2, self.nloc))
                arr[0,]=as_array(self.locx, 25602, shape=(self.nloc,))
                arr[1,]=as_array(self.locy, 25602, shape=(self.nloc,))
                return arr
    def free(self):
        lib.cellfree_do(byref(self))
class csc(Structure):#CSC sparse matrix
    _fields_=[
        ('id', c_uint32),
        ('x', c_void_p),
        ('nx', c_long),
        ('ny', c_long),
        ('header', c_char_p),
        ('nzmax', c_long),
        ('p', c_void_p),
        ('i', c_void_p),
        ('nref', c_void_p),
    ]
    def __init__(self, arr=None): #convert from numpy to C. Memory is borrowed
        dtype2id={#Conversion from sparse type to maos id
            np.float32:  25607,
            np.float64:  25601,
            np.complex64: 25606,
            np.complex128:25600,
        }
        if arr is not None and sp.isspmatrix_csc(arr):
            self.id=dtype2id.get(arr.dtype.type)
            #save subarrays
            self.xp=arr.data
            self.ip=arr.indices.astype(np.long)
            self.pp=arr.indptr.astype(np.long) #p
            self.x=self.xp.ctypes.data_as(c_void_p) #data
            self.i=self.ip.ctypes.data_as(c_void_p) #row index
            self.p=self.pp.ctypes.data_as(c_void_p)
            self.nx, self.ny=arr.shape #Fortran order
            self.nzmax=self.pp[-1]
        else:
            self.id=dtype2id.get(np.float64)
            self.x=None
            self.i=None
            self.p=None
            self.nx=0
            self.ny=0
            self.nzmax=0
        self.header=None
        self.nref=None
    def as_array(self): #convert form C to numpy. Memory is copied
        if self.nzmax>0:
            self.xp=as_array(self.x, self.id, (self.nzmax,))
            self.ip=as_array(self.i, 25603, (self.nzmax,))
            self.pp=as_array(self.p, 25603, (self.ny+1,))
            return sp.csc_matrix((self.xp, self.ip, self.pp), shape=(self.nx, self.ny))
        else:
            return sp.csc_matrix((self.nx,self.ny))
    def free(self):
        lib.cellfree_do(byref(self))
def convert_fields(fields):
    val2type={
        '*':c_void_p,
        'double':c_double,
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

#Create a ctypes class with field listed 
def make_class(name, fields):
    newfields=convert_fields(fields)
    class newclass(Structure):
        pass
        def as_array(self):#convert struct into dictionary
            out=dict()
            for ff in self._fields_:
                #convert C pointers to POINTER then to array
                if ff[1] is c_void_p:
                    exec('out[\''+ff[0]+'\']=cast(self.'+ff[0]+',POINTER(cell)).contents.as_array()')
                else:
                    exec('out[\''+ff[0]+'\']=self.'+ff[0])
            return out
        def free(self):
            print('to implement: free');
    newclass._fields_=newfields
    return newclass
