#!/usr/bin/env python3
"""Generate scripts/libaos.py, the ctypes Python bindings of aolib.so.

Usage:
    lib2py.py [srcdir [output]]

    srcdir: MAOS source tree (the directory containing maos/, lib/, mex/).
        Defaults to the parent directory of this script.
    output: file to write. Defaults to srcdir/scripts/libaos.py.

Inputs:
    Function declarations and struct definitions are parsed (with
    maos_parse.py) from srcdir/lib/*.h, srcdir/sys/scheduler_client.h and
    srcdir/mex/aolib.h. maos/parms.h and maos/types.h are not used.

Output:
    A Python module that does "from interface import *" (scripts/interface.py
    provides lib, the ctypes types, py2cell, ct2py, make_class, ...) and
    defines one wrapper per parsed C function that has arguments. For each
    function, the C prototype is used as the docstring; restype is set; numpy
    and Python values are converted to C types; pointer-to-pointer arguments
    (type**) and leading real*/double* arguments are treated as outputs and
    returned (after the C return value, if any) instead of taken as inputs.
    Structs used by the wrappers are recreated with make_class(). Functions
    with an argument or return type that cannot be mapped are skipped and
    reported on stdout.

The build does not run this script; run it by hand after changing the
headers (scripts/libaos.py is not tracked by git). Progress (one line per
function, plus unknown types and skipped functions) is printed to stdout.
"""
#2019-07-18 First version, derived from lib2mex.py

#from __future__ import print_function
import os, sys
import maos_parse
from pathlib import Path
import json
import glob
fnout=None
if len(sys.argv)>3 or (len(sys.argv)>1 and sys.argv[1] in ('-h','--help')):
    print('Usage: '+sys.argv[0]+' [srcdir [output]]', file=sys.stderr)
    sys.exit(0 if len(sys.argv)<=3 else 1)
if len(sys.argv)>1:
    srcdir=sys.argv[1]
    if len(sys.argv)>2:
        fnout=sys.argv[2]
else:
    srcdir=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if fnout is None:
    fnout=srcdir+'/scripts/libaos.py'
if not os.path.isdir(srcdir+'/maos/'):
    sys.exit('lib2py.py: '+srcdir+' is not the maos source dir (no maos/ subdirectory)')
if not os.path.isfile(srcdir+'/mex/aolib.h'):
    sys.exit('lib2py.py: '+srcdir+'/mex/aolib.h not found')
simu_all=list();
headerlist=glob.glob(srcdir+'/lib/*.h')
headerlist.append(srcdir+'/sys/scheduler_client.h')
headerlist.append(srcdir+'/mex/aolib.h')
structs=maos_parse.parse_structs('', headerlist)
funcs=maos_parse.parse_func('',headerlist)

verbose=1 #verbose message
#map between MAOS C types and python types
aotype2py={
    'dmat':('cell',0x6402),
    'lmat':('cell',0x6403),
    'cmat':('cell',0x6404),
    'smat':('cell',0x6408),
    'zmat':('cell',0x6409),
    'cell':('cell',0x6421),
    'map_t':('cell',0x6402),
    'loc_t':('loc',0x6402),
    'pts_t':('loc',0x6402),
    'dsp':('cell',0x6401),
    'csp':('cell',0x6400),
    'ssp':('cell',0x6407),
    'zsp':('cell',0x6406),
    'char':('c_char',0),
    'double':('c_double',0),
    'float':('c_float',0),
    'real':('c_double',0),
    'size_t':('c_size_t',0),
    'long':('c_long',0),
    'int':('c_int',0),
    'uint32_t':('c_uint',0),
    'uint64_t':('c_ulonglong',0),
    'void':('',0),
}
#process variable type
def handle_type(argtype, argname, leading):
    """Map one C function argument to Python wrapper code.

    Parameters:
        argtype: C type as returned by maos_parse.parse_func (e.g. "dmat*").
        argname: argument name, possibly with a "[n]" array suffix.
        leading: true if all preceding arguments were outputs. A real* or
            double* argument is then taken as a scalar output.

    Returns:
        (py2c, prep_type, prep_in, prep_out, pyarg_out):
        py2c: expression passed to the C function.
        prep_type: struct type name if the argument is a struct (its
            make_class definition must be emitted), else None or ''.
        prep_in: code line(s) converting the Python input before the call.
        prep_out: code line(s) allocating storage for an output argument.
            Non-empty means the argument is an output, not a Python input.
        pyarg_out: expression returning the output value to Python, or ''.
        If the type is not supported, every element (or py2c/argname) is
        'Unknown', which makes the caller skip the function.
    """
    if argname.count('[')>0 and argname[-1]==']' and (argtype=='real' or argtype=='double'):
        nel=argname[argname.find('[')+1:-1]
        if nel.isdigit(): #fixed size array is output
            argname=argname[:argname.find('[')]
            return (argname, None, '', f'    {argname}=(c_double*{nel})()\n', f'np.array({argname})')
    if argtype=='anyarray':
        argtype='cell*'
    elif argtype=='panyarray':
        argtype='cell**'

    if (len(argtype)>3 and argtype[-3:]=='***') or argname.count('[')>0: #unknown
        return ('Unknown', 'Unknown','Unknown','Unknown', 'Unknown')
    elif argtype[-2:]=='**': #output
        isout=1
        isref=1
        argtype=argtype[:-2]
    elif (argtype=='real*' or  argtype=='double*') and leading: #leading scalar pointers are output
        isout=2
        isref=1
        argtype=argtype[:-1]
    elif argtype[-1:]=='*': #pointer
        isout=0
        isref=1
        argtype=argtype[:-1]
    else:
        isout=0
        isref=0

    if argtype[-4:]=='cell':
        argtype='cell'

    pytype,tid=aotype2py.get(argtype, (None,None))
    prep_type=None #make classes
    prep_in=''  #convert input
    prep_out='' #allocation storage for C output 
    pyarg_out='' #return output value
    if pytype is None: #represent as a struct
        if structs.get(argtype,None) :
            #prep_type=argtype+'=make_class(\''+argtype+'\''+','+json.dumps(structs.get(argtype,None))+')'
            prep_type=argtype
            #py2c='pointer(dict_to_struct('+argname+', struct='+argtype+'))'
            py2c=argname
            prep_in=f'    if isinstance({argname}, dict): {argname}=pointer(dict_to_struct({argname}, struct={argtype})) \n'
        else:
            print('Unknown argument:', argtype, argname)
            pytype='Unknown'
            py2c='Unknown'
            argname='Unknown'
        
    elif pytype=='cell':#arbitrary array
        if isref:
            py2c='py2cellref({},{})'.format(argname, tid)
        else:
            py2c='py2cell({},{})'.format(argname, tid)
    else:
        py2c=pytype+'('+argname+')'
        if isref: #pointer input
            if pytype=='loc':
                py2c='byref('+py2c+')'
            elif (pytype=='c_double' or pytype=='c_float') and isout:
                py2c='byref('+argname+')'
            elif pytype=='c_double' or pytype=='c_float': #array input
                py2c='np.ascontiguousarray({},dtype={}).ctypes.data_as(c_void_p)'.format(argname, pytype)
            elif pytype=='c_char':
                py2c='c_char_p('+argname+'.encode("utf-8"))'
            else:
                py2c=argname+'.ctypes.data_as(c_void_p)'

    if isout==1: #pointer to pointer argument is an output (C allocates the object)
        prep_out=f'    {argname}=POINTER({pytype})()\n' #create memory for pointer to pointer (output)
        if prep_type:
            pyarg_out=argname #generic struct stay in ctypes
        else:
            pyarg_out='ct2py('+argname+')'
    elif isout==2: #function pointer input is output for scalar value
        prep_out=f'    {argname}={pytype}()\n'
        pyarg_out=argname+'.value'

    return (py2c, prep_type, prep_in, prep_out, pyarg_out)

#process function return type
def handle_output(funtype, funname):
    """Map a C function return type to Python wrapper code.

    Parameters:
        funtype: C return type (not 'void').
        funname: C function name.

    Returns:
        (fun_rtn, fun_val, str_def):
        fun_rtn: code line setting lib.<funname>.restype ('Unknown' in it if
            the type is not supported).
        fun_val: expression converting the result "out" to Python.
        str_def: struct type name if the result is a struct, else None.
    """
    if funtype[-1:]=='*': #pointer
        ispointer=1;
        funtype=funtype[:-1]
    else:
        ispointer=0
    if funtype[-4:]=='cell':
        funtype='cell'
    fun_rtn,fun_tid=aotype2py.get(funtype, (None,None))
    str_def=None
    if fun_rtn is None:#represent as a struct
        if structs.get(funtype, None):
            str_def=funtype
            #strdef[funtype]=structs[funtype]
            fun_rtn=funtype #'make_class(\''+funtype+'\''+','+json.dumps(structs.get(funtype,None))+')'
        else:
            print('Unknown output:', funtype, funname)
            fun_rtn='Unknown'
    if ispointer:
        if str_def:
            fun_val='out'
        else:
            fun_val='ct2py(out)'
        fun_rtn='POINTER('+fun_rtn+')'
    else:
        fun_val='out'
    if len(fun_rtn)==0:
        fun_rtn='None'
    return ('    lib.'+funname+'.restype='+fun_rtn, fun_val, str_def)

#Generate one Python wrapper per C function.
strdef={} #definition of structs
funcalls=list()
for funname in funcs: #loop over functions
    funtype=funcs[funname][0]  #function return type
    funargs=funcs[funname][1]  #Function arguments
    funname2=funcs[funname][2] #C function name
    if len(funargs)==0 :
        continue
    
    fundef='' #definition of functions

    #Define Python function
    pyargin='' #python function definition arguments
    argin=''   #arguments passing to C function
    pyargout=''  #output from python function
    prepin=''
    if funtype!='void': #with C outupt
        fun_rtn, fun_val, str_def=handle_output(funtype, funname2)
        pyargout+=fun_val+','
        argout='out='
        if structs.get(str_def, None):
            strdef[str_def]=structs[str_def]
    else:
        argout=''

    leading=1 #only output arguments so far
    for arg in funargs: #loop over the arguments
        argtype=arg[0]
        argname=arg[1]

        #Rename arguments that would shadow Python keywords, builtins, the
        #loc class of interface.py, or the local variable "out".
        basename=argname.split('[')[0]
        if basename in ['loc','input','in','out','str','string','len','format','type','min','max']:
            argname=basename+'_'+argname[len(basename):]
        if argtype=='number': #literal numbers
            argin+=argname+','
            leading=0
        else: #variables
            py2c, prep_type, prep_in, prep_out, pyarg_out=handle_type(argtype, argname, leading)
            if not prep_out:
                leading=0
            argin+=py2c+',' #input argument to C function
            if prep_type and structs.get(prep_type, None):
                strdef[prep_type]=structs.get(prep_type, None)
            if prep_in: #C function input preparation
                prepin+=prep_in
            if prep_out:
                prepin+=prep_out
            if pyarg_out:
                pyargout+=pyarg_out+','
            if not prep_out or len(funargs)==1: #when only 1 argument, assume is also input for Python function
                pyargin+=argname+','#python function input argument
    if verbose:
        print(funname, pyargin)
    if len(pyargin)>0 and pyargin[-1]==',':
        pyargin=pyargin[0:-1]
    if argin[-1]==',':
        argin=argin[0:-1]
    if len(pyargout)>0 and pyargout[-1]==',':
        pyargout=pyargout[0:-1]
    if len(prepin)>0 and prepin[-1]=='\n':
        prepin=prepin[0:-1]
    fundef+='def '+funname+'('+pyargin+'):'   #def function
    #C definition as documentation
    fundef+='\n    """'+funtype+' '+funname2+'('
    for arg in funargs:
        fundef+=arg[0]+' '+arg[1]+', '
    if len(funargs)>0:
        fundef=fundef[:-2]
    fundef+=');"""'
    if funtype!='void':
        fundef+='\n'+fun_rtn          #C function return type
    if len(prepin)>0:
        fundef+='\n'+prepin          #C function arguments
    fundef+='\n    '+argout+'lib.'+funname2+'('+argin+')'
    if len(pyargout)>0:
        fundef+='\n    return '+pyargout

    if(fundef.find('Unknown'))>-1:
        print('skip', funname)
        #print(fundef)
        #print("'''", file=fpout)
    else:
        funcalls.append(fundef)
        #print(fundef, file=fpout)
    #if(fundef.find('Unknown'))>-1:
    #    print("'''", file=fpout)

with open(fnout,'w') as fpout:
    print("#!/usr/bin/env python", file=fpout)
    print("\n#Do not modify! This file is automatically generated by lib2py.py.\n", file=fpout)
    print("from interface import *", file=fpout)
    for argtype, argval in strdef.items():
        print(argtype+'=make_class(\''+argtype+'\''+','+json.dumps(argval)+')', file=fpout)
    for fundef in funcalls:
        print(fundef, file=fpout)
